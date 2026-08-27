#!/usr/bin/env python3
"""Align assembled contigs and call haploid structural variants with SVIM-asm."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "svim_asm.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(f"Missing {CONFIG_FILE}; copy svim_asm.example.json to svim_asm.local.json.")
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    value = Path(config[key]).expanduser()
    return value if value.is_absolute() else PROJECT_ROOT / value


def sample_id(name):
    match = re.search(r"(\d{3}[A-Za-z])", name)
    return match.group(1).upper() if match else name.split("_")[0]


def run(command, stdout=None):
    print("Running:", " ".join(map(str, command)), flush=True)
    subprocess.run([str(value) for value in command], check=True, stdout=stdout)


def main():
    config = load_config()
    assemblies = path_value(config, "assemblies_dir")
    reference = path_value(config, "reference_fasta")
    output_root = path_value(config, "output_dir") / config.get("reference_name", "t2t")
    if not assemblies.is_dir():
        raise FileNotFoundError(f"Assembly directory does not exist: {assemblies}")
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
    output_root.mkdir(parents=True, exist_ok=True)
    work_root = output_root / "work"
    work_root.mkdir(exist_ok=True)
    suffix = config.get("assembly_suffix", "_p_ctg.fa")

    for assembly in sorted(assemblies.glob(f"*{suffix}")):
        sample = sample_id(assembly.name)
        sample_dir = output_root / sample
        final_vcf = sample_dir / "variants.vcf"
        if final_vcf.is_file() and final_vcf.stat().st_size > 0:
            print(f"Skipping completed sample: {sample}")
            continue
        sam = work_root / f"{sample}.sam"
        sorted_bam = work_root / f"{sample}.sorted.bam"
        with sam.open("wb") as output:
            run([
                config.get("minimap2", "minimap2"), "-a", "-x", "asm5", "--cs", "-r2k",
                "-t", config.get("threads", 32), reference, assembly,
            ], stdout=output)
        run([
            config.get("samtools", "samtools"), "sort", f"-m{config.get('sort_memory', '4G')}",
            f"-@{config.get('sort_threads', 4)}", "-o", sorted_bam, sam,
        ])
        run([config.get("samtools", "samtools"), "index", sorted_bam])
        run([config.get("svim_asm", "svim-asm"), "haploid", sample_dir, sorted_bam, reference])
        if final_vcf.is_file():
            sam.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
