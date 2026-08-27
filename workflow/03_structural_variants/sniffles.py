#!/usr/bin/env python3
"""Call structural variants independently for each aligned HiFi BAM with Sniffles."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "sniffles.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(f"Missing {CONFIG_FILE}; copy sniffles.example.json to sniffles.local.json.")
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    value = Path(config[key]).expanduser()
    return value if value.is_absolute() else PROJECT_ROOT / value


def sample_id(name):
    match = re.search(r"(\d{3}[A-Za-z])(?:\D|$)", name)
    return match.group(1).upper() if match else Path(name).stem


def main():
    config = load_config()
    bams = path_value(config, "aligned_bams")
    reference = path_value(config, "reference_fasta")
    repeats = path_value(config, "tandem_repeats")
    output_dir = path_value(config, "output_dir") / config.get("reference_name", "t2t")
    if not bams.is_dir():
        raise FileNotFoundError(f"Aligned BAM directory does not exist: {bams}")
    for label, path in (("reference FASTA", reference), ("tandem-repeat BED", repeats)):
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label}: {path}")
    output_dir.mkdir(parents=True, exist_ok=True)
    excluded = config.get("exclude_name_patterns", [])

    for bam in sorted(bams.glob("*.bam")):
        if any(pattern in bam.name for pattern in excluded):
            print(f"Skipping excluded input: {bam.name}")
            continue
        sample = sample_id(bam.name)
        output_vcf = output_dir / f"{sample}.vcf"
        if output_vcf.is_file() and output_vcf.stat().st_size > 0:
            print(f"Skipping completed sample: {sample}")
            continue
        command = [
            config.get("sniffles", "sniffles"), "--input", bam, "--vcf", output_vcf,
            "--reference", reference, "--allow-overwrite", "-t", config.get("threads", 32),
            "--tandem-repeats", repeats,
        ]
        print("Running:", " ".join(map(str, command)), flush=True)
        subprocess.run([str(value) for value in command], check=True)


if __name__ == "__main__":
    main()
