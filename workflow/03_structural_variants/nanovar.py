#!/usr/bin/env python3
"""Call structural variants from HiFi FASTQ files with NanoVar."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "nanovar.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(f"Missing {CONFIG_FILE}; copy nanovar.example.json to nanovar.local.json.")
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    value = Path(config[key]).expanduser()
    return value if value.is_absolute() else PROJECT_ROOT / value


def sample_id(name):
    match = re.search(r"(\d{3}[A-Za-z])(?:\D|$)", name)
    return match.group(1).upper() if match else name.split(".")[0]


def main():
    config = load_config()
    reads_dir = path_value(config, "hifi_reads")
    reference = path_value(config, "reference_fasta")
    result_subdir = config.get("result_subdir", config.get("reference_name", "reference"))
    output_root = path_value(config, "output_dir") / result_subdir
    if not reads_dir.is_dir():
        raise FileNotFoundError(f"HiFi FASTQ directory does not exist: {reads_dir}")
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
    output_root.mkdir(parents=True, exist_ok=True)
    suffix = config.get("fastq_suffix", ".fastq.gz")

    for reads in sorted(reads_dir.glob(f"*{suffix}")):
        sample = sample_id(reads.name)
        sample_dir = output_root / sample
        completed = list(sample_dir.glob("*.nanovar.pass.vcf")) if sample_dir.is_dir() else []
        if any(path.stat().st_size > 0 for path in completed):
            print(f"Skipping completed sample: {sample}")
            continue
        command = [
            config.get("nanovar", "nanovar"), reads, reference, sample_dir,
            "-t", config.get("threads", 32),
            "-l", config.get("minimum_sv_length", 50),
            "-x", config.get("read_type", "pacbio-ccs"),
        ]
        annotation = config.get("annotate_insertions")
        if annotation:
            command.extend(["--annotate_ins", annotation])
        print("Running:", " ".join(map(str, command)), flush=True)
        subprocess.run([str(value) for value in command], check=True)


if __name__ == "__main__":
    main()
