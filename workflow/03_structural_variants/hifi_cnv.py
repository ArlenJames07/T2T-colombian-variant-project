#!/usr/bin/env python3
"""Call copy-number variants from phased HiFi BAM and small-variant files."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "hifi_cnv.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(f"Missing {CONFIG_FILE}; copy hifi_cnv.example.json to hifi_cnv.local.json.")
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
    bams = path_value(config, "phased_bams")
    small_variants = path_value(config, "small_variants")
    reference = path_value(config, "reference_fasta")
    excluded_regions = path_value(config, "excluded_regions")
    output_root = path_value(config, "output_dir") / config.get("reference_name", "reference")
    for label, directory in (("phased BAM", bams), ("small-variant", small_variants)):
        if not directory.is_dir():
            raise FileNotFoundError(f"{label} directory does not exist: {directory}")
    for label, path in (("reference FASTA", reference), ("excluded-region BED", excluded_regions)):
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label}: {path}")
    output_root.mkdir(parents=True, exist_ok=True)

    bam_suffix = config.get("bam_suffix", ".bam")
    vcf_suffix = config.get("small_vcf_suffix", ".small.vcf.gz")
    for bam in sorted(bams.glob(f"*{bam_suffix}")):
        sample = sample_id(bam.name)
        sample_dir = output_root / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        prefix = sample_dir / sample
        final_vcf = Path(f"{prefix}.vcf.gz")
        if final_vcf.is_file() and final_vcf.stat().st_size > 0:
            print(f"Skipping completed sample: {sample}")
            continue
        maf_vcf = small_variants / f"{sample}{vcf_suffix}"
        if not maf_vcf.is_file():
            raise FileNotFoundError(f"Missing phased small-variant file: {maf_vcf}")
        command = [
            config.get("hificnv", "hificnv"), "--bam", bam, "--ref", reference,
            "--maf", maf_vcf, "--exclude", excluded_regions,
            "--threads", config.get("threads", 32), "--output-prefix", prefix,
        ]
        print("Running:", " ".join(map(str, command)), flush=True)
        subprocess.run([str(value) for value in command], check=True)


if __name__ == "__main__":
    main()
