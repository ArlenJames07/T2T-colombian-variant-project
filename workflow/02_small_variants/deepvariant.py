#!/usr/bin/env python3
"""Call PacBio HiFi SNVs and small indels with the DeepVariant container.

Aligned BAMs are discovered from the configured input directory. Outputs are
grouped by reference and named only with the canonical sample code, for example
``001P.vcf.gz`` and ``001P.g.vcf.gz``. Existing complete output sets are skipped.
"""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

try:
    import yaml
except ImportError as exc:  # pragma: no cover - depends on the runtime environment
    raise SystemExit("PyYAML is required: install it with 'python3 -m pip install pyyaml'") from exc


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BARCODE_PATTERN = re.compile(r"bc[0-9]+(?:v[0-9]+)?", re.IGNORECASE)
CODE_PATTERN = re.compile(r"_([0-9]{3}[A-Z])\.bam$", re.IGNORECASE)


def project_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file does not exist: {path}")
    with path.open(encoding="utf-8") as handle:
        content = yaml.safe_load(handle) or {}
    if not isinstance(content, dict):
        raise ValueError(f"Expected a YAML mapping in {path}")
    return content


def config_value(config: dict[str, Any], dotted_key: str, default: Any = None) -> Any:
    value: Any = config
    for key in dotted_key.split("."):
        if not isinstance(value, dict) or key not in value:
            if default is not None:
                return default
            raise KeyError(f"Missing required configuration value: {dotted_key}")
        value = value[key]
    return value


def completed(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def output_paths(output_dir: Path, code: str) -> tuple[Path, Path, Path, Path]:
    vcf = output_dir / f"{code}.vcf.gz"
    gvcf = output_dir / f"{code}.g.vcf.gz"
    return vcf, Path(f"{vcf}.tbi"), gvcf, Path(f"{gvcf}.tbi")


def display_command(command: Iterable[str | Path]) -> str:
    return shlex.join(str(part) for part in command)


def run(command: list[str | Path], *, dry_run: bool) -> None:
    print(f"[run] {display_command(command)}")
    if not dry_run:
        subprocess.run([str(part) for part in command], check=True)


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
            if not match:
                continue
            barcode = match.group(0).lower()
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paths", default="config/paths.yaml", help="Path configuration YAML.")
    parser.add_argument("--params", default="config/params.yaml", help="Parameter YAML.")
    parser.add_argument("--input", type=Path, help="Override inputs.aligned_bams.")
    parser.add_argument("--output", type=Path, help="Override small_variants.output root.")
    parser.add_argument("--reference", type=Path, help="Override references.fasta.")
    parser.add_argument("--reference-name", help="Reference-specific output parent name.")
    parser.add_argument("--threads", type=int, help="Override small_variants.threads.")
    parser.add_argument(
        "--sample", action="append", default=[], metavar="CODE",
        help="Only process sample codes containing CODE; repeat as needed.",
    )
    parser.add_argument(
        "--exclude", action="append", default=[], metavar="TEXT",
        help="Additionally exclude BAM filenames containing TEXT; repeat as needed.",
    )
    parser.add_argument("--force", action="store_true", help="Replace completed outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned work only.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        paths = load_yaml(project_path(args.paths))
        params = load_yaml(project_path(args.params))
        custom_reference = args.reference is not None
        reference = project_path(args.reference or config_value(paths, "references.fasta")).resolve()
        if not reference.is_file():
            raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
        configured_name = None if custom_reference else config_value(paths, "references.name", "t2t")
        parent_name = reference_name(reference, args.reference_name or configured_name)
        bam_dir = project_path(args.input or config_value(paths, "inputs.aligned_bams"))
        output_root = project_path(
            args.output or config_value(paths, "small_variants.output", "results/snvs")
        )
        output_dir = output_root / parent_name
        threads = args.threads or int(config_value(params, "small_variants.threads", 32))
        if threads < 1:
            raise ValueError("--threads must be at least 1")
        exclusions = [
            *config_value(params, "small_variants.exclude_name_patterns", []),
            *args.exclude,
        ]
        mapping = barcode_to_code(project_path(config_value(paths, "samples", "config/samples.tsv")))
        bams = find_bams(bam_dir, mapping, args.sample, exclusions)
        if not args.dry_run:
            output_dir.mkdir(parents=True, exist_ok=True)

        image = config_value(
            params, "small_variants.deepvariant_image", "google/deepvariant:1.8.0"
        )
        model_type = config_value(params, "small_variants.model_type", "PACBIO")
        ran = skipped = 0
        for code, bam in bams:
            vcf, vcf_index, gvcf, gvcf_index = output_paths(output_dir, code)
            expected = (vcf, vcf_index, gvcf, gvcf_index)
            if all(completed(path) for path in expected) and not args.force:
                print(f"[skip] {code}: complete DeepVariant outputs exist ({output_dir})")
                skipped += 1
                continue
            existing = [path for path in expected if path.exists() or path.is_symlink()]
            if existing and not args.force:
                names = ", ".join(path.name for path in existing)
                raise RuntimeError(
                    f"Partial DeepVariant output set for {code}: {names}. "
                    "Inspect the files, then use --force to replace the complete set."
                )
            if args.force and not args.dry_run:
                for path in expected:
                    path.unlink(missing_ok=True)
            command: list[str | Path] = [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{bam.parent}:/input:ro",
                "-v",
                f"{reference.parent}:/reference:ro",
                "-v",
                f"{output_dir}:/output",
                image,
                "/opt/deepvariant/bin/run_deepvariant",
                f"--model_type={model_type}",
                f"--ref=/reference/{reference.name}",
                f"--reads=/input/{bam.name}",
                f"--output_vcf=/output/{vcf.name}",
                f"--output_gvcf=/output/{gvcf.name}",
                f"--num_shards={threads}",
            ]
            run(command, dry_run=args.dry_run)
            if not args.dry_run and not all(completed(path) for path in expected):
                raise RuntimeError(f"DeepVariant did not produce all expected outputs for {code}")
            ran += 1
        action = "would run" if args.dry_run else "ran"
        print(f"[summary] DeepVariant {action}: {ran}; skipped complete samples: {skipped}")
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
