#!/usr/bin/env python3
"""Discover PacBio SV signatures, call SVs, and apply the legacy svpack filter."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "pbsv.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(f"Missing {CONFIG_FILE}; copy pbsv.example.json to pbsv.local.json.")
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    value = Path(config[key]).expanduser()
    return value if value.is_absolute() else PROJECT_ROOT / value


def sample_id(name):
    match = re.search(r"(\d{3}[A-Za-z])(?:\D|$)", name)
    return match.group(1).upper() if match else Path(name).stem


def complete(path):
    return path.is_file() and path.stat().st_size > 0


def run(command, stdout=None):
    print("Running:", " ".join(map(str, command)), flush=True)
    subprocess.run([str(value) for value in command], check=True, stdout=stdout)


def main():
    config = load_config()
    bams = path_value(config, "aligned_bams")
    reference = path_value(config, "reference_fasta")
    repeats = path_value(config, "tandem_repeats")
    root = path_value(config, "output_dir") / config.get("reference_name", "reference")
    if not bams.is_dir():
        raise FileNotFoundError(f"Aligned BAM directory does not exist: {bams}")
    for label, path in (("reference FASTA", reference), ("tandem-repeat BED", repeats)):
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label}: {path}")
    signatures = root / "signatures"
    raw_dir = root / "prefilter"
    filtered_dir = root / "filtered"
    for directory in (signatures, raw_dir, filtered_dir):
        directory.mkdir(parents=True, exist_ok=True)

    excluded = config.get("exclude_name_patterns", [])
    for bam in sorted(bams.glob("*.bam")):
        if any(pattern in bam.name for pattern in excluded):
            print(f"Skipping excluded input: {bam.name}")
            continue
        sample = sample_id(bam.name)
        signature = signatures / f"{sample}.svsig.gz"
        raw_vcf = raw_dir / f"{sample}.vcf"
        filtered_vcf = filtered_dir / f"{sample}.vcf"
        compressed_filtered = filtered_dir / f"{sample}.vcf.gz"

        if complete(filtered_vcf) or complete(compressed_filtered):
            print(f"Skipping completed sample: {sample}")
            continue
        if not complete(raw_vcf):
            if not complete(signature):
                run([config["pbsv"], "discover", bam, signature, "--tandem-repeats", repeats])
            else:
                print(f"Skipping signature discovery: {sample}")
            run([config["pbsv"], "call", reference, signature, "-j", config["threads"], raw_vcf])
        else:
            print(f"Skipping existing pbsv call: {sample}")

        temporary = filtered_vcf.with_suffix(".vcf.tmp")
        with temporary.open("wb") as output:
            run([
                config["svpack"], "filter", "--pass-only", "--min-svlen",
                config.get("minimum_sv_length", 50), raw_vcf,
            ], stdout=output)
        temporary.replace(filtered_vcf)


if __name__ == "__main__":
    main()
