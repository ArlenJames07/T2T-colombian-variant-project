#!/usr/bin/env python3
"""Call structural variants independently for each aligned HiFi BAM with cuteSV."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "cutesv.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(f"Missing {CONFIG_FILE}; copy cutesv.example.json to cutesv.local.json.")
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
    root = path_value(config, "output_dir") / config.get("reference_name", "t2t")
    if not bams.is_dir():
        raise FileNotFoundError(f"Aligned BAM directory does not exist: {bams}")
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
    raw_dir = root / "prefilter"
    work_dir = root / "work"
    raw_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    excluded = config.get("exclude_name_patterns", [])
    for bam in sorted(bams.glob("*.bam")):
        if any(pattern in bam.name for pattern in excluded):
            print(f"Skipping excluded input: {bam.name}")
            continue
        sample = sample_id(bam.name)
        output_vcf = raw_dir / f"{sample}.vcf"
        if output_vcf.is_file() and output_vcf.stat().st_size > 0:
            print(f"Skipping completed sample: {sample}")
            continue
        sample_work = work_dir / sample
        sample_work.mkdir(parents=True, exist_ok=True)
        command = [
            config.get("cutesv", "cuteSV"), bam, reference, output_vcf, sample_work,
            "--threads", config.get("threads", 32),
            "--min_support", config.get("minimum_support", 3),
            "--min_size", config.get("minimum_sv_length", 50),
            "--max_cluster_bias_INS", config.get("max_cluster_bias_ins", 1000),
            "--diff_ratio_merging_INS", config.get("diff_ratio_merging_ins", 0.9),
            "--max_cluster_bias_DEL", config.get("max_cluster_bias_del", 1000),
            "--diff_ratio_merging_DEL", config.get("diff_ratio_merging_del", 0.5),
        ]
        print("Running:", " ".join(map(str, command)), flush=True)
        subprocess.run([str(value) for value in command], check=True)


if __name__ == "__main__":
    main()
