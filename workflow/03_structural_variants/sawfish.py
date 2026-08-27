#!/usr/bin/env python3
"""Run Sawfish discovery/joint calling and optionally split SV and CNV records."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "sawfish.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(f"Missing {CONFIG_FILE}; copy sawfish.example.json to sawfish.local.json.")
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


def run(command):
    print("Running:", " ".join(map(str, command)), flush=True)
    subprocess.run([str(value) for value in command], check=True)


def split_joint_vcf(config, sample, joint_vcf, output_root):
    sv_dir = output_root / "SV_only"
    cnv_dir = output_root / "CNV_only"
    sv_dir.mkdir(exist_ok=True)
    cnv_dir.mkdir(exist_ok=True)
    sv_vcf = sv_dir / f"{sample}.SV_only.vcf.gz"
    cnv_vcf = cnv_dir / f"{sample}.CNV_depth_only.vcf.gz"

    if not complete(sv_vcf):
        run([config.get("bcftools", "bcftools"), "view", "-i", 'INFO/SVCLAIM~"J"',
             "-Oz", "-o", sv_vcf, joint_vcf])
        run([config.get("tabix", "tabix"), "-p", "vcf", sv_vcf])
    if not complete(cnv_vcf):
        run([config.get("bcftools", "bcftools"), "view", "-i", 'INFO/SVCLAIM="D"',
             "-Oz", "-o", cnv_vcf, joint_vcf])
        run([config.get("tabix", "tabix"), "-p", "vcf", cnv_vcf])


def main():
    config = load_config()
    bams = path_value(config, "phased_bams")
    reference = path_value(config, "reference_fasta")
    excluded_regions = path_value(config, "excluded_regions")
    output_root = path_value(config, "output_dir")
    if not bams.is_dir():
        raise FileNotFoundError(f"Phased BAM directory does not exist: {bams}")
    for label, path in (("reference FASTA", reference), ("excluded-region BED", excluded_regions)):
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label}: {path}")
    output_root.mkdir(parents=True, exist_ok=True)

    for bam in sorted(bams.glob("*.bam")):
        sample = sample_id(bam.name)
        discovery_dir = output_root / sample
        joint_dir = output_root / f"{sample}.joint"
        joint_vcf = joint_dir / "genotyped.sv.vcf.gz"
        if not complete(discovery_dir / "candidate.sv.bcf"):
            run([
                config.get("sawfish", "sawfish"), "discover", "--bam", bam,
                "--ref", reference, "--cnv-excluded-regions", excluded_regions,
                "--threads", config.get("threads", 32), "--output-dir", discovery_dir,
            ])
        else:
            print(f"Skipping Sawfish discovery: {sample}")
        if not complete(joint_vcf):
            run([
                config.get("sawfish", "sawfish"), "joint-call",
                "--threads", config.get("threads", 32), "--sample", discovery_dir,
                "--output-dir", joint_dir,
            ])
        else:
            print(f"Skipping Sawfish joint call: {sample}")
        if config.get("split_joint_vcfs", True) and complete(joint_vcf):
            split_joint_vcf(config, sample, joint_vcf, output_root)


if __name__ == "__main__":
    main()
