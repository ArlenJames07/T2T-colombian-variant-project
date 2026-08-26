#!/usr/bin/env python3
"""Scaffold hifiasm FASTAs with RagTag and extract chromosome-only FASTAs.

By default, this module scaffolds every ``<code>_<assembly-type>.fasta`` file in
``results/assemblies/fasta`` against the configured T2T reference. Results are
grouped by reference name and existing non-empty outputs are skipped.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, TextIO

try:
    import yaml
except ImportError as exc:  # pragma: no cover - depends on the runtime environment
    raise SystemExit("PyYAML is required: install it with 'python3 -m pip install pyyaml'") from exc


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ASSEMBLY_SUFFIXES = {
    "_hap1.fasta": "hap1",
    "_hap2.fasta": "hap2",
    "_p_ctg.fasta": "p_ctg",
}
DEFAULT_CHROMOSOME_REGEX = r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$"


def project_path(value: str | Path) -> Path:
    """Resolve a path relative to the repository, not the current directory."""
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


def display_command(command: Iterable[str | Path]) -> str:
    return shlex.join(str(part) for part in command)


def run(command: list[str | Path], *, dry_run: bool) -> None:
    print(f"[run] {display_command(command)}")
    if not dry_run:
        subprocess.run([str(part) for part in command], check=True)


def parse_assembly(path: Path) -> tuple[str, str] | None:
    for suffix, assembly_type in ASSEMBLY_SUFFIXES.items():
        if path.name.endswith(suffix):
            code = path.name[: -len(suffix)]
            return (code, assembly_type) if code else None
    return None


def find_assemblies(
    input_dir: Path,
    sample_filters: list[str],
    assembly_types: list[str],
) -> list[tuple[Path, str, str]]:
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Assembly input directory does not exist: {input_dir}")
    selected: list[tuple[Path, str, str]] = []
    seen: set[tuple[str, str]] = set()
    for path in sorted(input_dir.glob("*.fasta")):
        parsed = parse_assembly(path)
        if parsed is None:
            continue
        code, assembly_type = parsed
        if sample_filters and not any(term in code for term in sample_filters):
            continue
        if assembly_types and assembly_type not in assembly_types:
            continue
        key = (code, assembly_type)
        if key in seen:
            raise RuntimeError(f"Duplicate assembly code and type: {code}/{assembly_type}")
        seen.add(key)
        selected.append((path, code, assembly_type))
    if not selected:
        raise RuntimeError("No assembly FASTAs matched the requested filters")
    return selected


def reference_name(reference: Path, explicit_name: str | None) -> str:
    name = explicit_name
    if not name:
        name = reference.name
        for suffix in (".gz", ".bgz", ".fasta", ".fna", ".fa"):
            if name.lower().endswith(suffix):
                name = name[: -len(suffix)]
    if not name or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name):
        raise ValueError(
            "Reference name must start with a letter or number and contain only "
            "letters, numbers, dots, underscores, or hyphens"
        )
    return name


def open_fasta_text(path: Path) -> TextIO:
    if path.name.lower().endswith((".gz", ".bgz")):
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def reference_sequence_names(reference: Path) -> list[str]:
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
    names: list[str] = []
    with open_fasta_text(reference) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].strip().split(maxsplit=1)[0]
                if name:
                    names.append(name)
    if not names:
        raise ValueError(f"Reference contains no FASTA records: {reference}")
    return names


def load_chromosome_list(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Chromosome list does not exist: {path}")
    chromosomes: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                chromosomes.append(value.split()[0])
    if not chromosomes:
        raise ValueError(f"Chromosome list is empty: {path}")
    return chromosomes


def select_chromosomes(
    reference: Path,
    chromosome_list: Path | None,
    chromosome_regex: str,
) -> list[str]:
    reference_names = reference_sequence_names(reference)
    if chromosome_list:
        requested = load_chromosome_list(chromosome_list)
        unknown = sorted(set(requested).difference(reference_names))
        if unknown:
            raise ValueError(
                "Chromosome list contains IDs absent from the reference: " + ", ".join(unknown)
            )
        requested_set = set(requested)
        chromosomes = [name for name in reference_names if name in requested_set]
    else:
        try:
            pattern = re.compile(chromosome_regex)
        except re.error as exc:
            raise ValueError(f"Invalid --chromosome-regex: {exc}") from exc
        chromosomes = [name for name in reference_names if pattern.fullmatch(name)]
    if not chromosomes:
        raise ValueError(
            "No reference sequences were selected as chromosomes; provide "
            "--chromosome-regex or --chromosome-list"
        )
    return chromosomes


def normalized_ragtag_header(header: str, chromosomes: set[str]) -> str | None:
    if header in chromosomes:
        return header
    for suffix in ("_RagTag", "_ragtag"):
        if header.endswith(suffix) and header[: -len(suffix)] in chromosomes:
            return header[: -len(suffix)]
    return None


def extract_chromosomes(source: Path, destination: Path, chromosomes: list[str]) -> None:
    """Stream selected chromosome records to an atomic, normalized FASTA."""
    chromosome_set = set(chromosomes)
    found: set[str] = set()
    current: str | None = None
    written_bases = 0
    temporary = destination.with_suffix(f"{destination.suffix}.tmp")
    try:
        with source.open(encoding="utf-8") as input_handle, temporary.open(
            "w", encoding="utf-8"
        ) as output_handle:
            for line in input_handle:
                if line.startswith(">"):
                    header = line[1:].strip().split(maxsplit=1)[0]
                    current = normalized_ragtag_header(header, chromosome_set)
                    if current:
                        if current in found:
                            raise ValueError(
                                f"Scaffold FASTA contains chromosome {current} more than once: {source}"
                            )
                        found.add(current)
                        output_handle.write(f">{current}\n")
                elif current:
                    output_handle.write(line)
                    written_bases += len(line.strip())
        if not found or written_bases == 0:
            raise ValueError(f"No chromosome sequences were found in scaffold FASTA: {source}")
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise

    missing = [chromosome for chromosome in chromosomes if chromosome not in found]
    if missing:
        print(f"[warn] {destination.name}: chromosomes not present: {', '.join(missing)}")


def expose_scaffold(source: Path, destination: Path, *, force: bool) -> None:
    if not completed(source):
        raise RuntimeError(f"RagTag did not produce a non-empty scaffold FASTA: {source}")
    if destination.exists() or destination.is_symlink():
        if not force:
            if completed(destination):
                return
            raise FileExistsError(f"Refusing to replace existing result: {destination}")
        destination.unlink()
    relative_source = Path(os.path.relpath(source, destination.parent))
    destination.symlink_to(relative_source)


def scaffold_command(
    executable: Path,
    reference: Path,
    assembly: Path,
    output: Path,
    threads: int,
    parameters: dict[str, int | float],
    force: bool,
) -> list[str | Path]:
    command: list[str | Path] = [
        executable,
        "scaffold",
        "-t",
        str(threads),
        "-f",
        str(parameters["minimum_unique_alignment"]),
        "-q",
        str(parameters["minimum_mapq"]),
        "-d",
        str(parameters["maximum_merge_distance"]),
        "-i",
        str(parameters["minimum_grouping_confidence"]),
        "-a",
        str(parameters["minimum_location_confidence"]),
        "-s",
        str(parameters["minimum_orientation_confidence"]),
        "-o",
        output,
    ]
    if force:
        command.append("-w")
    command.extend([reference, assembly])
    return command


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paths", default="config/paths.yaml", help="Path configuration YAML.")
    parser.add_argument("--params", default="config/params.yaml", help="Parameter YAML.")
    parser.add_argument("--input", type=Path, help="Override assembly.fasta input directory.")
    parser.add_argument("--output", type=Path, help="Override assembly.scaffolds output root.")
    parser.add_argument("--reference", type=Path, help="Override references.fasta.")
    parser.add_argument(
        "--reference-name",
        help="Result parent directory name; defaults to configured name or reference stem.",
    )
    parser.add_argument("--threads", type=int, help="Override assembly.ragtag_threads.")
    parser.add_argument(
        "--sample", action="append", default=[], metavar="CODE",
        help="Only process sample codes containing CODE; repeat as needed.",
    )
    parser.add_argument(
        "--assembly-type",
        action="append",
        choices=("hap1", "hap2", "p_ctg"),
        default=[],
        help="Only process this assembly type; repeat as needed.",
    )
    parser.add_argument(
        "--chromosome-regex",
        help="Regex selecting reference chromosome IDs (default: human chr1-22, X and Y).",
    )
    parser.add_argument(
        "--chromosome-list",
        type=Path,
        help="File with one reference chromosome ID per line; overrides the regex.",
    )
    parser.add_argument("--min-unique-alignment", type=int)
    parser.add_argument("--min-mapq", type=int)
    parser.add_argument("--max-merge-distance", type=int)
    parser.add_argument("--min-grouping-confidence", type=float)
    parser.add_argument("--min-location-confidence", type=float)
    parser.add_argument("--min-orientation-confidence", type=float)
    parser.add_argument("--force", action="store_true", help="Replace completed outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned work only.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        paths = load_yaml(project_path(args.paths))
        params = load_yaml(project_path(args.params))
        custom_reference = args.reference is not None
        reference = project_path(args.reference or config_value(paths, "references.fasta"))
        configured_name = None if custom_reference else config_value(paths, "references.name", "t2t")
        parent_name = reference_name(reference, args.reference_name or configured_name)
        input_dir = project_path(args.input or config_value(paths, "assembly.fasta"))
        output_root = project_path(
            args.output or config_value(paths, "assembly.scaffolds", "results/assemblies/scaffolds")
        )
        threads = args.threads or int(config_value(params, "assembly.ragtag_threads", 32))
        if threads < 1:
            raise ValueError("--threads must be at least 1")

        chromosome_regex = args.chromosome_regex or config_value(
            params, "assembly.ragtag_chromosome_regex", DEFAULT_CHROMOSOME_REGEX
        )
        chromosome_list = project_path(args.chromosome_list) if args.chromosome_list else None
        chromosomes = select_chromosomes(reference, chromosome_list, chromosome_regex)
        print(f"[reference] {parent_name}: {reference}")
        print(f"[chromosomes] {len(chromosomes)} selected: {', '.join(chromosomes)}")

        ragtag_parameters: dict[str, int | float] = {
            "minimum_unique_alignment": args.min_unique_alignment
            if args.min_unique_alignment is not None
            else int(config_value(params, "assembly.ragtag.minimum_unique_alignment", 1000)),
            "minimum_mapq": args.min_mapq
            if args.min_mapq is not None
            else int(config_value(params, "assembly.ragtag.minimum_mapq", 10)),
            "maximum_merge_distance": args.max_merge_distance
            if args.max_merge_distance is not None
            else int(config_value(params, "assembly.ragtag.maximum_merge_distance", 100000)),
            "minimum_grouping_confidence": args.min_grouping_confidence
            if args.min_grouping_confidence is not None
            else float(config_value(params, "assembly.ragtag.minimum_grouping_confidence", 0.2)),
            "minimum_location_confidence": args.min_location_confidence
            if args.min_location_confidence is not None
            else float(config_value(params, "assembly.ragtag.minimum_location_confidence", 0.0)),
            "minimum_orientation_confidence": args.min_orientation_confidence
            if args.min_orientation_confidence is not None
            else float(config_value(params, "assembly.ragtag.minimum_orientation_confidence", 0.0)),
        }

        assemblies = find_assemblies(input_dir, args.sample, args.assembly_type)
        reference_root = output_root / parent_name
        full_dir = reference_root / "fasta"
        chromosome_dir = reference_root / "chromosomes"
        work_root = reference_root / "ragtag"
        if not args.dry_run:
            full_dir.mkdir(parents=True, exist_ok=True)
            chromosome_dir.mkdir(parents=True, exist_ok=True)
            work_root.mkdir(parents=True, exist_ok=True)

        executable = project_path(config_value(paths, "software.ragtag", "ragtag.py"))
        ragtag_ran = ragtag_skipped = filters_ran = filters_skipped = 0
        for assembly, code, assembly_type in assemblies:
            stem = f"{code}_{assembly_type}"
            work_dir = work_root / stem
            ragtag_scaffold = work_dir / "ragtag.scaffold.fasta"
            full_result = full_dir / f"{stem}.fasta"
            chromosome_result = chromosome_dir / f"{stem}.fasta"

            if completed(full_result) and not args.force:
                print(f"[skip] {stem}: scaffold FASTA exists ({full_result})")
                ragtag_skipped += 1
            elif completed(ragtag_scaffold) and not args.force:
                print(f"[link] {ragtag_scaffold} -> {full_result}")
                if not args.dry_run:
                    expose_scaffold(ragtag_scaffold, full_result, force=False)
                ragtag_skipped += 1
            else:
                if not args.dry_run:
                    work_dir.mkdir(parents=True, exist_ok=True)
                command = scaffold_command(
                    executable=executable,
                    reference=reference,
                    assembly=assembly,
                    output=work_dir,
                    threads=threads,
                    parameters=ragtag_parameters,
                    force=args.force,
                )
                run(command, dry_run=args.dry_run)
                if not args.dry_run:
                    expose_scaffold(ragtag_scaffold, full_result, force=args.force)
                ragtag_ran += 1

            if completed(chromosome_result) and not args.force:
                print(f"[skip] {stem}: chromosome FASTA exists ({chromosome_result})")
                filters_skipped += 1
            else:
                print(f"[filter] {full_result} -> {chromosome_result}")
                if not args.dry_run:
                    extract_chromosomes(full_result, chromosome_result, chromosomes)
                filters_ran += 1

        action = "would run" if args.dry_run else "ran"
        print(
            f"[summary] RagTag {action}: {ragtag_ran}; RagTag skipped: {ragtag_skipped}; "
            f"chromosome filters {action}: {filters_ran}; filters skipped: {filters_skipped}"
        )
        return 0
    except (
        FileExistsError,
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
