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


def pass_filter_complete(vcf):
    index = Path(f"{vcf}.tbi")
    marker = Path(f"{vcf}.pass_only.done")
    return (
        complete(vcf)
        and complete(index)
        and marker.is_file()
        and marker.stat().st_mtime_ns >= vcf.stat().st_mtime_ns
    )


def keep_pass_variants(config, vcf):
    """Atomically replace a VCF with its PASS records and rebuild its index."""
    temporary_vcf = vcf.with_name(f".{vcf.name}.pass.tmp.vcf.gz")
    temporary_index = Path(f"{temporary_vcf}.tbi")
    final_index = Path(f"{vcf}.tbi")
    marker = Path(f"{vcf}.pass_only.done")

    for stale_path in (temporary_vcf, temporary_index):
        stale_path.unlink(missing_ok=True)

    run([
        config.get("bcftools", "bcftools"), "view",
        "--apply-filters", "PASS",
        "--output-type", "z",
        "--output", temporary_vcf,
        vcf,
    ])
    run([
        config.get("tabix", "tabix"), "--force", "--preset", "vcf",
        temporary_vcf,
    ])
    temporary_vcf.replace(vcf)
    temporary_index.replace(final_index)
    marker.write_text("FILTER=PASS\n", encoding="utf-8")
    print(f"Kept PASS variants and removed RefCall records: {vcf}")


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

    excluded = config.get("exclude_name_patterns", [])
    for bam in sorted(aligned_bams.glob("*.bam")):
        if any(pattern in bam.name for pattern in excluded):
            print(f"Skipping excluded input: {bam.name}")
            continue
        sample = sample_id(bam.name)
        vcf = output_dir / f"{sample}.vcf.gz"
        gvcf = output_dir / f"{sample}.g.vcf.gz"
        if pass_filter_complete(vcf):
            print(f"Skipping completed sample: {sample}")
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

        keep_pass_variants(config, vcf)


if __name__ == "__main__":
    main()
