#!/usr/bin/env python3
"""Run or resume HiFi assembly and convert primary-contig GFA files to FASTA.

Paths and thread counts come from the repository configuration by default. Existing
non-empty outputs are skipped unless ``--force`` is used, so the workflow can be
resumed without repeating completed, expensive work.
"""

from __future__ import annotations

import argparse
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
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
RAW_GFA_SUFFIXES = {
    ".asm.bp.hap1.p_ctg.gfa": "hap1",
    ".asm.bp.hap2.p_ctg.gfa": "hap2",
    ".asm.bp.p_ctg.gfa": "p_ctg",
}
NORMALIZED_GFA_SUFFIXES = {
    "_hap1.gfa": "hap1",
    "_hap2.gfa": "hap2",
    "_p_ctg.gfa": "p_ctg",
}


def project_path(value: str | Path) -> Path:
    """Resolve a path relative to the repository, independently of the current directory."""
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


def strip_fastq_suffix(filename: str) -> str:
    for suffix in FASTQ_SUFFIXES:
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    raise ValueError(f"Unsupported FASTQ filename: {filename}")


def sample_name(read: Path) -> str:
    """Map 01_<sample>.fastq.gz to the established <sample> output name."""
    return strip_fastq_suffix(read.name).removeprefix("01_")


def sample_code(name: str) -> str:
    """Return the cohort code at the end of an established sample name."""
    code = name.rsplit("_", maxsplit=1)[-1]
    if not code:
        raise ValueError(f"Could not derive a sample code from: {name}")
    return code


def selected(name: str, include: list[str], exclude: list[str]) -> bool:
    return (not include or any(term in name for term in include)) and not any(
        term in name for term in exclude
    )


def completed(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def display_command(command: Iterable[str | Path]) -> str:
    return shlex.join(str(part) for part in command)


def run(command: list[str | Path], *, dry_run: bool, stdout: Any = None) -> None:
    print(f"[run] {display_command(command)}")
    if not dry_run:
        subprocess.run([str(part) for part in command], check=True, stdout=stdout)


def find_reads(reads_dir: Path) -> list[Path]:
    if not reads_dir.is_dir():
        raise FileNotFoundError(f"HiFi read directory does not exist: {reads_dir}")
    return sorted(
        path
        for path in reads_dir.iterdir()
        if path.is_file() and path.name.endswith(FASTQ_SUFFIXES)
    )


def assemble(
    reads_dir: Path,
    work_dir: Path,
    gfa_dir: Path,
    hifiasm: Path,
    threads: int,
    include: list[str],
    exclude: list[str],
    force: bool,
    dry_run: bool,
) -> tuple[int, int]:
    reads = [read for read in find_reads(reads_dir) if selected(read.name, include, exclude)]
    if not reads:
        raise RuntimeError("No FASTQ files matched the requested sample filters")
    samples: dict[str, tuple[str, Path]] = {}
    for read in reads:
        name = sample_name(read)
        code = sample_code(name)
        if code in samples:
            other = samples[code][1]
            raise RuntimeError(
                f"Multiple FASTQ files resolve to sample code {code}: {other} and {read}. "
                "Select one run with --sample or exclude one in config/params.yaml."
            )
        samples[code] = (name, read)

    if not dry_run:
        work_dir.mkdir(parents=True, exist_ok=True)
        gfa_dir.mkdir(parents=True, exist_ok=True)

    ran = skipped = 0
    for code, (name, read) in sorted(samples.items()):
        prefix = work_dir / f"02_{name}.asm"
        marker = gfa_dir / f"{code}_p_ctg.gfa"
        if completed(marker) and not force:
            print(f"[skip] {code}: complete GFA exists ({marker})")
            skipped += 1
            continue
        run([hifiasm, "-o", prefix, "-t", str(threads), read], dry_run=dry_run)
        if not dry_run:
            expose_primary_gfas(prefix, gfa_dir, code, force=force)
        ran += 1
    return ran, skipped


def expose_primary_gfas(prefix: Path, gfa_dir: Path, code: str, *, force: bool) -> None:
    """Expose hifiasm primary-contig outputs under stable, short result names."""
    for suffix, assembly_type in RAW_GFA_SUFFIXES.items():
        source = Path(f"{prefix}{suffix}")
        destination = gfa_dir / f"{code}_{assembly_type}.gfa"
        if not completed(source):
            raise RuntimeError(f"Expected hifiasm output was not produced: {source}")
        if destination.exists() or destination.is_symlink():
            if not force:
                raise FileExistsError(f"Refusing to replace existing GFA: {destination}")
            destination.unlink()
        destination.symlink_to(source.resolve())


def parse_gfa(path: Path) -> tuple[str, str] | None:
    for suffix, assembly_type in NORMALIZED_GFA_SUFFIXES.items():
        if path.name.endswith(suffix):
            return path.name[: -len(suffix)], assembly_type
    return None


def convert(
    gfa_dir: Path,
    fasta_dir: Path,
    gfatools: Path,
    include: list[str],
    exclude: list[str],
    force: bool,
    dry_run: bool,
) -> tuple[int, int]:
    if not gfa_dir.is_dir():
        raise FileNotFoundError(f"GFA directory does not exist: {gfa_dir}")
    gfas = [
        (path, parsed)
        for path in sorted(gfa_dir.glob("*.gfa"))
        if (parsed := parse_gfa(path)) is not None
        and selected(parsed[0], include, exclude)
    ]
    if not gfas:
        raise RuntimeError("No primary-contig GFA files matched the requested sample filters")
    if not dry_run:
        fasta_dir.mkdir(parents=True, exist_ok=True)

    ran = skipped = 0
    for gfa, (name, assembly_type) in gfas:
        fasta = fasta_dir / f"{name}_{assembly_type}.fasta"
        if completed(fasta) and not force:
            print(f"[skip] {name}/{assembly_type}: FASTA exists ({fasta})")
            skipped += 1
            continue

        command: list[str | Path] = [gfatools, "gfa2fa", gfa]
        print(f"[run] {display_command(command)} > {shlex.quote(str(fasta))}")
        if not dry_run:
            temporary = fasta.with_suffix(f"{fasta.suffix}.tmp")
            try:
                with temporary.open("wb") as output_handle:
                    subprocess.run(
                        [str(part) for part in command], check=True, stdout=output_handle
                    )
                if temporary.stat().st_size == 0:
                    raise RuntimeError(f"gfatools produced an empty FASTA: {temporary}")
                temporary.replace(fasta)
            except BaseException:
                temporary.unlink(missing_ok=True)
                raise
        ran += 1
    return ran, skipped


def executable_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else project_path(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--step", choices=("assemble", "convert", "all"), default="all",
        help="Pipeline stage to run (default: all).",
    )
    parser.add_argument("--paths", default="config/paths.yaml", help="Path configuration YAML.")
    parser.add_argument("--params", default="config/params.yaml", help="Parameter YAML.")
    parser.add_argument("--reads", type=Path, help="Override inputs.hifi_reads.")
    parser.add_argument("--work-output", type=Path, help="Override assembly.hifiasm_work.")
    parser.add_argument("--gfa-output", type=Path, help="Override assembly.hifiasm_gfa.")
    parser.add_argument("--fasta-output", type=Path, help="Override assembly.fasta.")
    parser.add_argument("--threads", type=int, help="Override assembly.hifiasm_threads.")
    parser.add_argument(
        "--sample", action="append", default=[], metavar="TEXT",
        help="Only process sample codes containing TEXT; repeat to select multiple samples.",
    )
    parser.add_argument(
        "--exclude", action="append", default=[], metavar="TEXT",
        help="Additionally exclude filenames containing TEXT; repeat as needed.",
    )
    parser.add_argument("--force", action="store_true", help="Replace completed outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned work without running tools.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        paths = load_yaml(project_path(args.paths))
        params = load_yaml(project_path(args.params))

        reads_dir = project_path(args.reads or config_value(paths, "inputs.hifi_reads"))
        work_dir = project_path(
            args.work_output
            or config_value(paths, "assembly.hifiasm_work", "work/assembly/hifiasm")
        )
        gfa_dir = project_path(args.gfa_output or config_value(paths, "assembly.hifiasm_gfa"))
        fasta_dir = project_path(args.fasta_output or config_value(paths, "assembly.fasta"))
        threads = args.threads or int(config_value(params, "assembly.hifiasm_threads", 32))
        configured_exclusions = config_value(params, "assembly.exclude_name_patterns", [])
        exclude = [*configured_exclusions, *args.exclude]
        if threads < 1:
            raise ValueError("--threads must be at least 1")

        total_ran = total_skipped = 0
        if args.step in {"assemble", "all"}:
            ran, skipped = assemble(
                reads_dir=reads_dir,
                work_dir=work_dir,
                gfa_dir=gfa_dir,
                hifiasm=executable_path(config_value(paths, "software.hifiasm")),
                threads=threads,
                include=args.sample,
                exclude=exclude,
                force=args.force,
                dry_run=args.dry_run,
            )
            total_ran += ran
            total_skipped += skipped
        if args.step in {"convert", "all"}:
            ran, skipped = convert(
                gfa_dir=gfa_dir,
                fasta_dir=fasta_dir,
                gfatools=executable_path(config_value(paths, "software.gfatools")),
                include=args.sample,
                exclude=exclude,
                force=args.force,
                dry_run=args.dry_run,
            )
            total_ran += ran
            total_skipped += skipped

        action = "would run" if args.dry_run else "ran"
        print(f"[summary] {action}: {total_ran}; skipped existing: {total_skipped}")
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
