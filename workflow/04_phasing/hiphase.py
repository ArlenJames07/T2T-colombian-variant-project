#!/usr/bin/env python3
"""Phase small and structural variants together with HiPhase."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "hiphase.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(f"Missing {CONFIG_FILE}; copy hiphase.example.json to hiphase.local.json.")
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    value = Path(config[key]).expanduser()
    return value if value.is_absolute() else PROJECT_ROOT / value


def sample_id(name):
    match = re.search(r"(\d{3}[A-Za-z])(?:\D|$)", name)
    return match.group(1).upper() if match else name.split(".")[0]


def file_index(directory, pattern):
    index = {}
    for path in sorted(directory.glob(pattern)):
        index.setdefault(sample_id(path.name), path)
    return index


def complete(path):
    return path.is_file() and path.stat().st_size > 0


def run(command):
    print("Running:", " ".join(map(str, command)), flush=True)
    subprocess.run([str(value) for value in command], check=True)


def prepare_small_vcf(config, source_vcf, compressed_dir):
    """Return a BGZF-compressed, indexed VCF suitable for HiPhase."""
    if source_vcf.name.endswith(".vcf.gz"):
        return source_vcf

    compressed_vcf = compressed_dir / f"{source_vcf.name}.gz"
    compressed_index = Path(f"{compressed_vcf}.tbi")
    if (
        complete(compressed_vcf)
        and complete(compressed_index)
        and compressed_vcf.stat().st_mtime_ns >= source_vcf.stat().st_mtime_ns
        and compressed_index.stat().st_mtime_ns >= compressed_vcf.stat().st_mtime_ns
    ):
        return compressed_vcf

    temporary_vcf = compressed_dir / f".{source_vcf.name}.tmp.vcf.gz"
    temporary_index = Path(f"{temporary_vcf}.tbi")
    for stale_path in (temporary_vcf, temporary_index):
        stale_path.unlink(missing_ok=True)

    try:
        run([
            config.get("bcftools", "bcftools"), "view", "--output-type", "z",
            "--output", temporary_vcf, source_vcf,
        ])
        run([
            config.get("tabix", "tabix"), "--force", "--preset", "vcf",
            temporary_vcf,
        ])
        temporary_vcf.replace(compressed_vcf)
        temporary_index.replace(compressed_index)
    except Exception:
        temporary_vcf.unlink(missing_ok=True)
        temporary_index.unlink(missing_ok=True)
        raise

    print(f"Prepared BGZF VCF for HiPhase: {compressed_vcf}")
    return compressed_vcf


def main():
    config = load_config()
    small_variants = path_value(config, "small_variants")
    bams = path_value(config, "aligned_bams")
    structural_variants = path_value(config, "structural_variants")
    reference = path_value(config, "reference_fasta")
    output_root = path_value(config, "output_dir")
    bam_output = output_root / "bamfiles"
    variant_output = output_root / "variants"
    compressed_input = output_root / "compressed_inputs"
    for label, directory in (
        ("small-variant", small_variants), ("aligned BAM", bams),
        ("structural-variant", structural_variants),
    ):
        if not directory.is_dir():
            raise FileNotFoundError(f"{label} directory does not exist: {directory}")
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
    bam_output.mkdir(parents=True, exist_ok=True)
    variant_output.mkdir(parents=True, exist_ok=True)
    compressed_input.mkdir(parents=True, exist_ok=True)

    bam_by_sample = file_index(bams, "*.bam")
    small_by_sample = {}
    for path in sorted(small_variants.glob(config.get("small_vcf_glob", "*.vcf.gz"))):
        if not path.name.endswith(".g.vcf.gz"):
            small_by_sample.setdefault(sample_id(path.name), path)
    sv_by_sample = file_index(structural_variants, f"*{config.get('sv_suffix', '.vcf.gz')}")

    for sample, sv_vcf in sv_by_sample.items():
        bam = bam_by_sample.get(sample)
        small_vcf = small_by_sample.get(sample)
        if bam is None or small_vcf is None:
            missing = "BAM" if bam is None else "small-variant VCF"
            print(f"Skipping {sample}: missing {missing}")
            continue

        phased_bam = bam_output / f"{sample}.bam"
        phased_small = variant_output / f"{sample}.small.vcf.gz"
        phased_sv = variant_output / f"{sample}.SV.vcf.gz"
        if all(complete(path) for path in (phased_bam, phased_small, phased_sv)):
            print(f"Skipping completed sample: {sample}")
            continue

        hiphase_small_vcf = prepare_small_vcf(config, small_vcf, compressed_input)
        command = [
            config["hiphase"], "--threads", config.get("threads", 32),
            "--reference", reference, "--bam", bam, "--output-bam", phased_bam,
            "--vcf", hiphase_small_vcf, "--output-vcf", phased_small,
            "--vcf", sv_vcf, "--output-vcf", phased_sv,
            "--stats-file", variant_output / f"{sample}.stats.csv",
            "--blocks-file", variant_output / f"{sample}.blocks.tsv",
            "--summary-file", variant_output / f"{sample}.summary.tsv",
        ]
        run(command)


if __name__ == "__main__":
    main()
