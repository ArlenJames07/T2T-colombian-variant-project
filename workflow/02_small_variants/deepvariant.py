#!/usr/bin/env python3
"""Call small variants from aligned HiFi reads with DeepVariant in Docker."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "deepvariant.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        example = CONFIG_FILE.with_name("deepvariant.example.json")
        raise FileNotFoundError(f"Missing {CONFIG_FILE}. Copy {example.name} and edit its paths.")
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    path = Path(config[key]).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def sample_id(filename):
    match = re.search(r"(\d{3}[A-Za-z])(?:\D|$)", filename)
    return match.group(1).upper() if match else Path(filename).stem


def complete(path):
    return path.is_file() and path.stat().st_size > 0


def run(command):
    print("Running:", " ".join(map(str, command)), flush=True)
    subprocess.run([str(value) for value in command], check=True)


def pass_filter_complete(source_vcf, filtered_vcf, compressed_filtered_vcf):
    """Return True when both PASS-only VCF outputs are current."""
    return (
        complete(source_vcf)
        and complete(filtered_vcf)
        and complete(compressed_filtered_vcf)
        and filtered_vcf.stat().st_mtime_ns >= source_vcf.stat().st_mtime_ns
        and compressed_filtered_vcf.stat().st_mtime_ns >= source_vcf.stat().st_mtime_ns
    )


def keep_pass_variants(config, source_vcf, filtered_vcf, compressed_filtered_vcf):
    """Write PASS records as plain and compressed VCFs without changing the source."""
    temporary_vcf = filtered_vcf.with_name(f".{filtered_vcf.name}.tmp")
    temporary_compressed_vcf = compressed_filtered_vcf.with_name(
        f".{compressed_filtered_vcf.name}.tmp"
    )
    for temporary in (temporary_vcf, temporary_compressed_vcf):
        temporary.unlink(missing_ok=True)

    try:
        run([
            config.get("bcftools", "bcftools"), "view",
            "--apply-filters", "PASS",
            "--output-type", "v",
            "--output", temporary_vcf,
            source_vcf,
        ])
        run([
            config.get("bcftools", "bcftools"), "view",
            "--output-type", "z",
            "--output", temporary_compressed_vcf,
            temporary_vcf,
        ])
        temporary_vcf.replace(filtered_vcf)
        temporary_compressed_vcf.replace(compressed_filtered_vcf)
    except Exception:
        for temporary in (temporary_vcf, temporary_compressed_vcf):
            temporary.unlink(missing_ok=True)
        raise

    print(f"Saved PASS-only variants: {filtered_vcf}, {compressed_filtered_vcf}")


def main():
    config = load_config()
    aligned_bams = path_value(config, "aligned_bams").resolve()
    reference = path_value(config, "reference_fasta").resolve()
    output_dir = (path_value(config, "output_dir") / config.get("reference_name", "reference")).resolve()
    if not aligned_bams.is_dir():
        raise FileNotFoundError(f"Aligned BAM directory does not exist: {aligned_bams}")
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
    output_dir.mkdir(parents=True, exist_ok=True)
    filtered_dir = output_dir / "filtered"
    filtered_dir.mkdir(exist_ok=True)

    excluded = config.get("exclude_name_patterns", [])
    for bam in sorted(aligned_bams.glob("*.bam")):
        if any(pattern in bam.name for pattern in excluded):
            print(f"Skipping excluded input: {bam.name}")
            continue
        sample = sample_id(bam.name)
        vcf = output_dir / f"{sample}.vcf.gz"
        gvcf = output_dir / f"{sample}.g.vcf.gz"
        filtered_vcf = filtered_dir / f"{sample}.vcf"
        compressed_filtered_vcf = filtered_dir / f"{sample}.vcf.gz"
        if pass_filter_complete(vcf, filtered_vcf, compressed_filtered_vcf):
            print(f"Skipping completed PASS-only VCF: {sample}")
            continue

        if not complete(vcf):
            run([
                config.get("docker", "docker"), "run", "--rm",
                "-v", f"{aligned_bams}:/input:ro",
                "-v", f"{reference.parent}:/reference:ro",
                "-v", f"{output_dir}:/output",
                config.get("image", "google/deepvariant:1.8.0"),
                "/opt/deepvariant/bin/run_deepvariant",
                f"--model_type={config.get('model_type', 'PACBIO')}",
                f"--ref=/reference/{reference.name}",
                f"--reads=/input/{bam.name}",
                f"--output_vcf=/output/{vcf.name}",
                f"--output_gvcf=/output/{gvcf.name}",
                f"--num_shards={config.get('threads', 32)}",
            ])

        keep_pass_variants(config, vcf, filtered_vcf, compressed_filtered_vcf)


if __name__ == "__main__":
    main()
