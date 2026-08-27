#!/usr/bin/env python3
"""Run or resume per-sample cuteSV calling and publish filtered SV VCFs.

Aligned BAMs are mapped to canonical sample codes through ``config/samples.tsv``.
The raw caller output is published under ``prefilter/``. PASS variants with
``|SVLEN| >= minimum_size`` when SVLEN is defined (plus BND calls without
SVLEN) are sorted, compressed, indexed, and published under ``filtered/``.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import config_value, load_yaml, project_path, run


BARCODE_PATTERN = re.compile(r"bc[0-9]+(?:v[0-9]+)?", re.IGNORECASE)
CODE_PATTERN = re.compile(r"_([0-9]{3}[A-Z])\.bam$", re.IGNORECASE)


def completed(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def reference_name(reference: Path, explicit_name: str | None) -> str:
    name = explicit_name
    if not name:
        name = reference.name
        for suffix in (".gz", ".bgz", ".fasta", ".fna", ".fa"):
            if name.lower().endswith(suffix):
                name = name[: -len(suffix)]
    if not name or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name):
        raise ValueError("Invalid reference name for an output directory")
    return name


def barcode_to_code(sample_table: Path) -> dict[str, str]:
    if not sample_table.is_file():
        raise FileNotFoundError(f"Sample table does not exist: {sample_table}")
    mapping: dict[str, str] = {}
    with sample_table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "bam"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Sample table must contain columns {sorted(required)}: {sample_table}")
        for row in reader:
            barcode_match = BARCODE_PATTERN.search(row["bam"])
            if not barcode_match:
                continue
            barcode = barcode_match.group(0).lower()
            code = row["sample_id"].strip()
            if barcode in mapping and mapping[barcode] != code:
                raise ValueError(f"Barcode {barcode} maps to multiple sample codes")
            mapping[barcode] = code
    return mapping


def canonical_code(bam: Path, mapping: dict[str, str]) -> str:
    barcode_match = BARCODE_PATTERN.search(bam.name)
    if barcode_match:
        mapped = mapping.get(barcode_match.group(0).lower())
        if mapped:
            return mapped
    code_match = CODE_PATTERN.search(bam.name)
    if not code_match:
        raise ValueError(f"Could not derive a sample code from BAM filename: {bam.name}")
    return code_match.group(1).upper()


def find_bams(
    bam_dir: Path,
    mapping: dict[str, str],
    sample_filters: list[str],
    exclusions: list[str],
) -> list[tuple[str, Path]]:
    if not bam_dir.is_dir():
        raise FileNotFoundError(f"Aligned BAM directory does not exist: {bam_dir}")
    selected: list[tuple[str, Path]] = []
    seen: dict[str, Path] = {}
    for bam in sorted(bam_dir.glob("*.bam")):
        if any(term in bam.name for term in exclusions):
            continue
        code = canonical_code(bam, mapping)
        if sample_filters and not any(term in code for term in sample_filters):
            continue
        if code in seen:
            raise RuntimeError(
                f"Multiple BAMs resolve to sample code {code}: {seen[code]} and {bam}. "
                "Select or exclude one sequencing run."
            )
        seen[code] = bam
        selected.append((code, bam))
    if not selected:
        raise RuntimeError("No aligned BAMs matched the requested filters")
    return selected


def output_paths(output_dir: Path, sample: str) -> tuple[Path, Path, Path]:
    prefilter_vcf = output_dir / "prefilter" / f"{sample}.vcf"
    filtered_vcf = output_dir / "filtered" / f"{sample}.vcf.gz"
    return prefilter_vcf, filtered_vcf, Path(f"{filtered_vcf}.tbi")


def prepare_for_replacement(path: Path, *, force: bool, dry_run: bool) -> None:
    if force and not dry_run and path.exists():
        path.unlink()


def call_sample(
    *,
    sample: str,
    bam: Path,
    reference: Path,
    output_dir: Path,
    work_dir: Path,
    cutesv: str,
    svpack: str,
    bcftools: str,
    tabix: str,
    threads: int,
    minimum_support: int,
    minimum_size: int,
    force: bool,
    dry_run: bool,
) -> bool:
    prefilter_vcf, filtered_vcf, filtered_index = output_paths(output_dir, sample)
    temporary_filtered_vcf = work_dir / "filtered" / f"{sample}.vcf"
    caller_work_dir = work_dir / "caller" / sample
    published = (prefilter_vcf, filtered_vcf, filtered_index)
    if not force and all(completed(path) for path in published):
        print(f"[skip] {sample}: all published cuteSV stages exist ({output_dir})")
        return False

    for path in (*published, temporary_filtered_vcf):
        prepare_for_replacement(path, force=force, dry_run=dry_run)
        if not dry_run:
            path.parent.mkdir(parents=True, exist_ok=True)
    if not dry_run:
        caller_work_dir.mkdir(parents=True, exist_ok=True)

    changed = False
    if completed(prefilter_vcf):
        print(f"[skip] {sample}: prefilter VCF exists ({prefilter_vcf})")
    else:
        run(
            [
                cutesv,
                bam,
                reference,
                prefilter_vcf,
                caller_work_dir,
                "--sample",
                sample,
                "--threads",
                threads,
                "--min_support",
                minimum_support,
                "--min_size",
                minimum_size,
                "--max_cluster_bias_INS",
                1000,
                "--diff_ratio_merging_INS",
                0.9,
                "--max_cluster_bias_DEL",
                1000,
                "--diff_ratio_merging_DEL",
                0.5,
            ],
            dry_run=dry_run,
        )
        changed = True

    if completed(filtered_vcf):
        print(f"[skip] {sample}: filtered VCF exists ({filtered_vcf})")
    else:
        if completed(temporary_filtered_vcf):
            print(f"[skip] {sample}: temporary filtered VCF exists ({temporary_filtered_vcf})")
        else:
            run(
                [svpack, "filter", "--pass-only", "--min-svlen", minimum_size, prefilter_vcf],
                stdout=temporary_filtered_vcf,
                dry_run=dry_run,
            )
        run(
            [bcftools, "sort", "-Oz", "-o", filtered_vcf, temporary_filtered_vcf],
            dry_run=dry_run,
        )
        changed = True

    if completed(filtered_index):
        print(f"[skip] {sample}: filtered VCF index exists ({filtered_index})")
    else:
        run([tabix, "-f", "-p", "vcf", filtered_vcf], dry_run=dry_run)
        changed = True
    return changed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "bam",
        nargs="?",
        type=Path,
        help="Optional single aligned BAM; otherwise discover BAMs from the configured input.",
    )
    parser.add_argument("--paths", default="config/paths.yaml", help="Path configuration YAML.")
    parser.add_argument("--params", default="config/params.yaml", help="Parameter YAML.")
    parser.add_argument("--input", type=Path, help="Override inputs.aligned_bams.")
    parser.add_argument("--output", type=Path, help="Override the published cuteSV stage root.")
    parser.add_argument("--work-output", type=Path, help="Override the cuteSV work root.")
    parser.add_argument("--reference", type=Path, help="Override references.fasta.")
    parser.add_argument("--reference-name", help="Reference-specific output parent name.")
    parser.add_argument("--threads", type=int, help="Override structural_variants.threads.")
    parser.add_argument("--min-support", type=int, help="Minimum supporting reads.")
    parser.add_argument("--min-size", type=int, help="Minimum absolute SVLEN (inclusive).")
    parser.add_argument(
        "--sample",
        action="append",
        default=[],
        metavar="CODE",
        help="Only process sample codes containing CODE; repeat as needed.",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        metavar="TEXT",
        help="Additionally exclude BAM filenames containing TEXT; repeat as needed.",
    )
    parser.add_argument("--force", action="store_true", help="Replace completed outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned work only.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.bam and args.input:
            raise ValueError("Use either the positional BAM or --input, not both")
        paths = load_yaml(project_path(args.paths))
        params = load_yaml(project_path(args.params))
        custom_reference = args.reference is not None
        reference = project_path(args.reference or config_value(paths, "references.fasta")).resolve()
        if not reference.is_file():
            raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
        configured_name = None if custom_reference else config_value(paths, "references.name", "t2t")
        parent_name = reference_name(reference, args.reference_name or configured_name)
        output_root = project_path(
            args.output
            or config_value(paths, "structural_variants.cutesv_output", "results/sv/cutesv")
        )
        work_root = project_path(
            args.work_output
            or config_value(paths, "structural_variants.cutesv_work", "work/structural_variants/cutesv")
        )
        output_dir = output_root / parent_name
        work_dir = work_root / parent_name
        threads = (
            args.threads
            if args.threads is not None
            else int(config_value(params, "structural_variants.threads", 32))
        )
        minimum_support = (
            args.min_support
            if args.min_support is not None
            else int(config_value(params, "structural_variants.cute_sv_min_support", 3))
        )
        minimum_size = (
            args.min_size
            if args.min_size is not None
            else int(config_value(params, "structural_variants.minimum_size", 50))
        )
        if threads < 1:
            raise ValueError("--threads must be at least 1")
        if minimum_support < 1:
            raise ValueError("--min-support must be at least 1")
        if minimum_size < 1:
            raise ValueError("--min-size must be at least 1")

        mapping = barcode_to_code(project_path(config_value(paths, "samples", "config/samples.tsv")))
        exclusions = [
            *config_value(params, "structural_variants.exclude_name_patterns", []),
            *config_value(params, "structural_variants.cutesv_exclude_name_patterns", []),
            *args.exclude,
        ]
        if args.bam:
            bam = project_path(args.bam).resolve()
            if not bam.is_file():
                raise FileNotFoundError(f"Aligned BAM does not exist: {bam}")
            code = canonical_code(bam, mapping)
            if args.sample and not any(term in code for term in args.sample):
                raise RuntimeError(f"BAM resolves to {code}, which does not match --sample")
            selected_bams = [(code, bam)]
        else:
            bam_dir = project_path(args.input or config_value(paths, "inputs.aligned_bams"))
            selected_bams = find_bams(bam_dir, mapping, args.sample, exclusions)

        cutesv = str(config_value(paths, "software.cutesv", "cuteSV"))
        svpack = str(config_value(paths, "software.svpack", "svpack"))
        bcftools = str(config_value(paths, "software.bcftools", "bcftools"))
        tabix = str(config_value(paths, "software.tabix", "tabix"))
        ran = skipped = 0
        for sample, bam in selected_bams:
            changed = call_sample(
                sample=sample,
                bam=bam,
                reference=reference,
                output_dir=output_dir,
                work_dir=work_dir,
                cutesv=cutesv,
                svpack=svpack,
                bcftools=bcftools,
                tabix=tabix,
                threads=threads,
                minimum_support=minimum_support,
                minimum_size=minimum_size,
                force=args.force,
                dry_run=args.dry_run,
            )
            if changed:
                ran += 1
            else:
                skipped += 1
        action = "would process" if args.dry_run else "processed"
        print(f"[summary] cuteSV {action}: {ran}; skipped complete samples: {skipped}")
        return 0
    except (
        FileNotFoundError,
        KeyError,
        RuntimeError,
        ValueError,
        subprocess.CalledProcessError,
    ) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
