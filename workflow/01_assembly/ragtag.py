#!/usr/bin/env python3
"""Scaffold each assembly against the configured reference with RagTag."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "ragtag.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        example = CONFIG_FILE.with_name("ragtag.example.json")
        raise FileNotFoundError(f"Missing {CONFIG_FILE}. Copy {example.name} and edit its paths.")
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    path = Path(config[key]).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def complete(path):
    return path.is_file() and path.stat().st_size > 0


def assembly_name(path, aliases):
    name = path.name
    for old, new in aliases.items():
        name = name.replace(old, new)
    for suffix in (".fasta", ".fa", ".fna"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    sample = re.search(r"(\d{3}[A-Za-z])", name)
    kind = next((kind for kind in ("hap1", "hap2", "p_ctg") if kind in name), None)
    return f"{sample.group(1).upper()}_{kind}" if sample and kind else name


def main():
    config = load_config()
    assemblies = path_value(config, "assemblies_dir")
    reference = path_value(config, "reference_fasta")
    output_root = path_value(config, "output_dir") / config.get("reference_name", "reference")
    ragtag_root = output_root / "ragtag"
    fasta_root = output_root / "fasta"
    if not assemblies.is_dir():
        raise FileNotFoundError(f"Assembly directory does not exist: {assemblies}")
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
    ragtag_root.mkdir(parents=True, exist_ok=True)
    fasta_root.mkdir(parents=True, exist_ok=True)

    suffix = config.get("fasta_suffix", ".fasta")
    for assembly in sorted(assemblies.glob(f"*{suffix}")):
        name = assembly_name(assembly, config.get("sample_aliases", {}))
        sample_dir = ragtag_root / name
        scaffold = sample_dir / "ragtag.scaffold.fasta"
        if complete(scaffold):
            print(f"Skipping scaffolded assembly: {name}")
        else:
            sample_dir.mkdir(parents=True, exist_ok=True)
            command = [
                config["ragtag"], "scaffold", reference, assembly,
                "-o", sample_dir, "-t", config["threads"],
            ]
            print("Running:", " ".join(map(str, command)), flush=True)
            subprocess.run([str(value) for value in command], check=True)

        stable_link = fasta_root / f"{name}.fasta"
        if complete(scaffold) and not stable_link.exists() and not stable_link.is_symlink():
            stable_link.symlink_to(Path("..") / "ragtag" / name / scaffold.name)


if __name__ == "__main__":
    main()
