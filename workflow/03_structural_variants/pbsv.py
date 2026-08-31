#!/usr/bin/env python3
"""Discover PacBio SV signatures, call SVs, and apply the legacy svpack filter."""

import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "pbsv.local.json"
MAX_CONCURRENT_FILES = 4


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


def filtered_outputs_complete(raw_vcf, filtered_vcf, compressed_vcf, index):
    """Return True when all filtered outputs required by HiPhase are current."""
    return (
        complete(raw_vcf)
        and complete(filtered_vcf)
        and complete(compressed_vcf)
        and complete(index)
        and filtered_vcf.stat().st_mtime_ns >= raw_vcf.stat().st_mtime_ns
        and compressed_vcf.stat().st_mtime_ns >= filtered_vcf.stat().st_mtime_ns
        and index.stat().st_mtime_ns >= compressed_vcf.stat().st_mtime_ns
    )


def filter_and_prepare_vcf(config, raw_vcf, filtered_vcf, compressed_vcf):
    """Filter SVs and create plain, BGZF-compressed, and indexed VCF outputs."""
    index = Path(f"{compressed_vcf}.tbi")
    temporary_vcf = filtered_vcf.with_name(f".{filtered_vcf.stem}.tmp.vcf")
    temporary_compressed = compressed_vcf.with_name(
        f".{filtered_vcf.stem}.tmp.vcf.gz"
    )
    temporary_index = Path(f"{temporary_compressed}.tbi")
    for temporary in (temporary_vcf, temporary_compressed, temporary_index):
        temporary.unlink(missing_ok=True)

    try:
        with temporary_vcf.open("wb") as output:
            run([
                config["svpack"], "filter", "--pass-only", "--min-svlen",
                config.get("minimum_sv_length", 50), raw_vcf,
            ], stdout=output)
        run([
            config.get("bcftools", "bcftools"), "view", "--output-type", "z",
            "--output", temporary_compressed, temporary_vcf,
        ])
        run([
            config.get("tabix", "tabix"), "--force", "--preset", "vcf",
            temporary_compressed,
        ])
        temporary_vcf.replace(filtered_vcf)
        temporary_compressed.replace(compressed_vcf)
        temporary_index.replace(index)
    except Exception:
        for temporary in (temporary_vcf, temporary_compressed, temporary_index):
            temporary.unlink(missing_ok=True)
        raise

    print(f"Prepared HiPhase inputs: {filtered_vcf}, {compressed_vcf}, {index}")


def process_bam(config, bam, reference, repeats, signatures, raw_dir, filtered_dir):
    """Run the complete pbsv workflow for one BAM file."""
    sample = sample_id(bam.name)
    signature = signatures / f"{sample}.svsig.gz"
    raw_vcf = raw_dir / f"{sample}.vcf"
    filtered_vcf = filtered_dir / f"{sample}.vcf"
    compressed_filtered = filtered_dir / f"{sample}.vcf.gz"
    compressed_index = Path(f"{compressed_filtered}.tbi")

    if filtered_outputs_complete(
        raw_vcf, filtered_vcf, compressed_filtered, compressed_index
    ):
        print(f"Skipping completed sample: {sample}")
        return
    if not complete(raw_vcf):
        if not complete(signature):
            run([
                config["pbsv"], "discover", bam, signature,
                "--tandem-repeats", repeats,
            ])
        else:
            print(f"Skipping signature discovery: {sample}")
        run([
            config["pbsv"], "call", reference, signature,
            "-j", config["threads"], raw_vcf,
        ])
    else:
        print(f"Skipping existing pbsv call: {sample}")

    filter_and_prepare_vcf(config, raw_vcf, filtered_vcf, compressed_filtered)


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
    input_bams = []
    for bam in sorted(bams.glob("*.bam")):
        if any(pattern in bam.name for pattern in excluded):
            print(f"Skipping excluded input: {bam.name}")
        else:
            input_bams.append(bam)

    print(f"Processing up to {MAX_CONCURRENT_FILES} BAM files concurrently")
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENT_FILES) as executor:
        futures = [
            executor.submit(
                process_bam,
                config,
                bam,
                reference,
                repeats,
                signatures,
                raw_dir,
                filtered_dir,
            )
            for bam in input_bams
        ]
        for future in futures:
            future.result()


if __name__ == "__main__":
    main()
