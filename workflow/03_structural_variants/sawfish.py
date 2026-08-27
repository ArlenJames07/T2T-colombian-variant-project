#!/usr/bin/env python3
"""Run, resume, or import Sawfish using its native per-sample directory layout.

Each sample is published as ``results/sv/sawfish/<code>/`` with its matching
joint-call directory at ``<code>.joint/``. A separate ``bedgraphs/<code>/``
folder collects the discovery and joint copy-number bedGraph tracks. Completed
Sawfish directories can be reused with ``--import-root`` without recomputation.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import config_value, load_yaml, project_path, run


BARCODE_PATTERN = re.compile(r"bc[0-9]+(?:v[0-9]+)?", re.IGNORECASE)
CODE_PATTERN = re.compile(r"_([0-9]{3}[A-Z])(?:\.bam)?$", re.IGNORECASE)


@dataclass(frozen=True)
class SampleJob:
    code: str
    bam: Path | None
    discovery: Path
    joint: Path
    import_discovery: Path | None = None
    import_joint: Path | None = None


def completed(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


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
            match = BARCODE_PATTERN.search(row["bam"])
            if match:
                mapping[match.group(0).lower()] = row["sample_id"].strip()
    return mapping


def canonical_code(path: Path, mapping: dict[str, str]) -> str:
    barcode = BARCODE_PATTERN.search(path.name)
    if barcode and barcode.group(0).lower() in mapping:
        return mapping[barcode.group(0).lower()]
    name = path.name[:-4] if path.name.endswith(".bam") else path.name
    match = CODE_PATTERN.search(name)
    if not match:
        raise ValueError(f"Could not derive a canonical sample code from: {path.name}")
    return match.group(1).upper()


def is_selected(code: str, name: str, samples: list[str], exclusions: list[str]) -> bool:
    return not any(term in name for term in exclusions) and (
        not samples or any(term in code for term in samples)
    )


def import_jobs(
    import_root: Path,
    output_root: Path,
    mapping: dict[str, str],
    samples: list[str],
    exclusions: list[str],
) -> list[SampleJob]:
    if not import_root.is_dir():
        raise FileNotFoundError(f"Sawfish import root does not exist: {import_root}")
    jobs: list[SampleJob] = []
    seen: set[str] = set()
    for source_discovery in sorted(import_root.iterdir()):
        if not source_discovery.is_dir() or source_discovery.name.endswith(".joint"):
            continue
        if not completed(source_discovery / "candidate.sv.bcf"):
            continue
        source_joint = import_root / f"{source_discovery.name}.joint"
        if not completed(source_joint / "genotyped.sv.vcf.gz"):
            raise RuntimeError(f"Missing completed joint folder for {source_discovery}")
        code = canonical_code(source_discovery, mapping)
        if not is_selected(code, source_discovery.name, samples, exclusions):
            continue
        if code in seen:
            raise RuntimeError(f"Multiple completed Sawfish runs resolve to {code}")
        seen.add(code)
        jobs.append(
            SampleJob(
                code=code,
                bam=None,
                discovery=output_root / code,
                joint=output_root / f"{code}.joint",
                import_discovery=source_discovery.resolve(),
                import_joint=source_joint.resolve(),
            )
        )
    if not jobs:
        raise RuntimeError("No completed Sawfish samples matched the requested filters")
    return jobs


def bam_jobs(
    bam: Path | None,
    bam_root: Path,
    output_root: Path,
    mapping: dict[str, str],
    samples: list[str],
    exclusions: list[str],
) -> list[SampleJob]:
    if bam is None and not bam_root.is_dir():
        raise FileNotFoundError(f"Aligned BAM directory does not exist: {bam_root}")
    alignments = [bam.resolve()] if bam else sorted(bam_root.glob("*.bam"))
    jobs: list[SampleJob] = []
    seen: set[str] = set()
    for alignment in alignments:
        if not alignment.is_file():
            raise FileNotFoundError(f"Aligned BAM does not exist: {alignment}")
        code = canonical_code(alignment, mapping)
        if not is_selected(code, alignment.name, samples, exclusions):
            continue
        if code in seen:
            raise RuntimeError(f"Multiple BAMs resolve to {code}")
        seen.add(code)
        jobs.append(
            SampleJob(
                code=code,
                bam=alignment,
                discovery=output_root / code,
                joint=output_root / f"{code}.joint",
            )
        )
    if not jobs:
        raise RuntimeError("No aligned BAMs matched the requested filters")
    return jobs


def publish_directory(source: Path, target: Path, *, force: bool, dry_run: bool) -> bool:
    if target.is_symlink() and target.resolve() == source.resolve():
        return False
    if target.is_symlink() and force:
        if not dry_run:
            target.unlink()
    elif target.exists() or target.is_symlink():
        raise RuntimeError(f"Cannot publish {source}: target already exists: {target}")
    print(f"$ ln -s {source} {target}")
    if not dry_run:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(source, target_is_directory=True)
    return True


def publish_file(source: Path, target: Path, *, force: bool, dry_run: bool) -> bool:
    if not source.is_file():
        raise FileNotFoundError(f"Expected Sawfish bedGraph does not exist: {source}")
    if target.is_symlink() and target.resolve() == source.resolve():
        return False
    if target.is_symlink() and force:
        if not dry_run:
            target.unlink()
    elif target.exists() or target.is_symlink():
        raise RuntimeError(f"Cannot publish {source}: target already exists: {target}")
    print(f"$ ln -s {source} {target}")
    if not dry_run:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(source)
    return True


def publish_bedgraphs(job: SampleJob, output_root: Path, *, force: bool, dry_run: bool) -> bool:
    sample_root = output_root / "bedgraphs" / job.code
    discovery = job.import_discovery or job.discovery
    joint = job.import_joint or job.joint
    changed = publish_file(
        discovery / "copynum.bedgraph",
        sample_root / "discovery.copynum.bedgraph",
        force=force,
        dry_run=dry_run,
    )
    sample_dirs = sorted((joint / "samples").glob("sample*"))
    if len(sample_dirs) != 1:
        raise RuntimeError(f"Expected one Sawfish sample directory under {joint}")
    changed |= publish_file(
        sample_dirs[0] / "copynum.bedgraph",
        sample_root / "joint.copynum.bedgraph",
        force=force,
        dry_run=dry_run,
    )
    return changed


def process_sample(
    job: SampleJob,
    *,
    reference: Path,
    cnv_excluded: Path,
    output_root: Path,
    sawfish: str,
    threads: int,
    force: bool,
    dry_run: bool,
) -> bool:
    changed = False
    if job.import_discovery and job.import_joint:
        changed |= publish_directory(
            job.import_discovery, job.discovery, force=force, dry_run=dry_run
        )
        changed |= publish_directory(job.import_joint, job.joint, force=force, dry_run=dry_run)
    else:
        if force or not completed(job.discovery / "candidate.sv.bcf"):
            if job.bam is None:
                raise RuntimeError(f"No BAM supplied for {job.code}")
            run(
                [
                    sawfish, "discover", "--bam", job.bam, "--ref", reference,
                    "--cnv-excluded-regions", cnv_excluded, "--threads", threads,
                    "--output-dir", job.discovery,
                ],
                dry_run=dry_run,
            )
            changed = True
        if force or not completed(job.joint / "genotyped.sv.vcf.gz"):
            run(
                [
                    sawfish, "joint-call", "--threads", threads, "--sample", job.discovery,
                    "--output-dir", job.joint,
                ],
                dry_run=dry_run,
            )
            changed = True
    changed |= publish_bedgraphs(
        job, output_root, force=force, dry_run=dry_run
    )
    if not changed:
        print(f"[skip] {job.code}: native Sawfish folders and bedGraphs exist")
    return changed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam", nargs="?", type=Path, help="Optional single aligned BAM.")
    parser.add_argument("--paths", default="config/paths.yaml", help="Path configuration YAML.")
    parser.add_argument("--params", default="config/params.yaml", help="Parameter YAML.")
    parser.add_argument("--input", type=Path, help="Override inputs.aligned_bams.")
    parser.add_argument("--import-root", type=Path, help="Reuse completed Sawfish folders.")
    parser.add_argument("--output", type=Path, help="Override results/sv/sawfish.")
    parser.add_argument("--reference", type=Path, help="Override references.fasta.")
    parser.add_argument("--threads", type=int, help="Override structural_variants.threads.")
    parser.add_argument("--sample", action="append", default=[], metavar="CODE")
    parser.add_argument("--exclude", action="append", default=[], metavar="TEXT")
    parser.add_argument("--force", action="store_true", help="Replace published links.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned work only.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.import_root and (args.bam or args.input):
            raise ValueError("--import-root cannot be combined with a BAM or --input")
        paths = load_yaml(project_path(args.paths))
        params = load_yaml(project_path(args.params))
        output_root = project_path(
            args.output or config_value(paths, "structural_variants.sawfish_output", "results/sv/sawfish")
        )
        reference = project_path(args.reference or config_value(paths, "references.fasta")).resolve()
        cnv_excluded = project_path(config_value(paths, "references.cnv_excluded")).resolve()
        if not reference.is_file():
            raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
        if not cnv_excluded.is_file():
            raise FileNotFoundError(f"CNV exclusion BED does not exist: {cnv_excluded}")
        threads = (
            args.threads
            if args.threads is not None
            else int(config_value(params, "structural_variants.threads", 32))
        )
        if threads < 1:
            raise ValueError("--threads must be at least 1")
        mapping = barcode_to_code(project_path(config_value(paths, "samples", "config/samples.tsv")))
        exclusions = [
            *config_value(params, "structural_variants.exclude_name_patterns", []),
            *args.exclude,
        ]
        if args.import_root:
            jobs = import_jobs(
                project_path(args.import_root), output_root, mapping, args.sample, exclusions
            )
        else:
            bam_root = project_path(args.input or config_value(paths, "inputs.aligned_bams"))
            jobs = bam_jobs(args.bam, bam_root, output_root, mapping, args.sample, exclusions)
        sawfish = str(config_value(paths, "software.sawfish", "sawfish"))
        changed = skipped = 0
        for job in jobs:
            if process_sample(
                job, reference=reference, cnv_excluded=cnv_excluded,
                output_root=output_root, sawfish=sawfish, threads=threads,
                force=args.force, dry_run=args.dry_run,
            ):
                changed += 1
            else:
                skipped += 1
        action = "would process" if args.dry_run else "processed"
        print(f"[summary] Sawfish {action}: {changed}; skipped complete samples: {skipped}")
        return 0
    except (FileNotFoundError, KeyError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
