#!/usr/bin/env python3
"""Build HiFi assemblies and convert the main hifiasm graphs to FASTA."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "hifiasm.local.json"
GFA_TYPES = {
    "bp.p_ctg.gfa": "p_ctg",
    "bp.hap1.p_ctg.gfa": "hap1",
    "bp.hap2.p_ctg.gfa": "hap2",
}


def load_config():
    if not CONFIG_FILE.is_file():
        example = CONFIG_FILE.with_name("hifiasm.example.json")
        raise FileNotFoundError(f"Missing {CONFIG_FILE}. Copy {example.name} and edit its paths.")
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    path = Path(config[key]).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def complete(path):
    return path.is_file() and path.stat().st_size > 0


def sample_id(name, aliases):
    for old, new in aliases.items():
        name = name.replace(old, new)
    match = re.search(r"(\d{3}[A-Za-z])(?:\D|$)", name)
    return match.group(1).upper() if match else name.split(".")[0]


def run(command, stdout=None):
    print("Running:", " ".join(map(str, command)), flush=True)
    subprocess.run([str(value) for value in command], check=True, stdout=stdout)


def index_bams(config, raw_bams):
    for bam in sorted(raw_bams.glob("*.bam")):
        index = Path(f"{bam}.pbi")
        if complete(index):
            print(f"Skipping indexed BAM: {bam.name}")
            continue
        run([config["pbindex"], bam])


def convert_bams(config, raw_bams, fastq_dir):
    fastq_dir.mkdir(parents=True, exist_ok=True)
    for bam in sorted(raw_bams.glob("*.bam")):
        prefix = fastq_dir / f"01_{bam.stem}"
        fastq = Path(f"{prefix}.fastq.gz")
        if complete(fastq):
            print(f"Skipping existing FASTQ: {fastq.name}")
            continue
        run([config["bam2fastq"], "-o", prefix, bam, "-j", config["threads"]])


def canonical_gfas(gfa_dir, sample):
    return {kind: gfa_dir / f"{sample}_{kind}.gfa" for kind in GFA_TYPES.values()}


def assemble_reads(config, fastq_dir, gfa_dir):
    gfa_dir.mkdir(parents=True, exist_ok=True)
    aliases = config.get("sample_aliases", {})
    excluded = config.get("exclude_name_patterns", [])
    suffix = config.get("fastq_suffix", ".fastq.gz")
    for fastq in sorted(fastq_dir.glob(f"*{suffix}")):
        if any(pattern in fastq.name for pattern in excluded):
            print(f"Skipping excluded input: {fastq.name}")
            continue
        sample = sample_id(fastq.name, aliases)
        outputs = canonical_gfas(gfa_dir, sample)
        if all(complete(path) for path in outputs.values()):
            print(f"Skipping assembled sample: {sample}")
            continue

        prefix = gfa_dir / f"{sample}.asm"
        run([config["hifiasm"], "-o", prefix, "-t", config["threads"], fastq])
        for suffix_name, kind in GFA_TYPES.items():
            generated = Path(f"{prefix}.{suffix_name}")
            destination = outputs[kind]
            if generated.exists() and not complete(destination):
                generated.replace(destination)


def convert_gfas(config, gfa_dir, fasta_dir):
    fasta_dir.mkdir(parents=True, exist_ok=True)
    for gfa in sorted(gfa_dir.glob("*.gfa")):
        if not any(gfa.name.endswith(f"_{kind}.gfa") for kind in GFA_TYPES.values()):
            continue
        fasta = fasta_dir / f"{gfa.stem}.fasta"
        if complete(fasta):
            print(f"Skipping existing FASTA: {fasta.name}")
            continue
        temporary = fasta.with_suffix(f"{fasta.suffix}.tmp")
        with temporary.open("wb") as output:
            run([config["gfatools"], "gfa2fa", gfa], stdout=output)
        temporary.replace(fasta)


def main():
    config = load_config()
    raw_bams = path_value(config, "raw_bams")
    fastq_dir = path_value(config, "fastq_dir")
    output_dir = path_value(config, "output_dir")
    gfa_dir = output_dir / "gfa"
    fasta_dir = output_dir / "fasta"

    if (config.get("index_bams", False) or config.get("convert_bams", False)) and not raw_bams.is_dir():
        raise FileNotFoundError(f"Raw BAM directory does not exist: {raw_bams}")
    if config.get("assemble", True) and not fastq_dir.is_dir():
        raise FileNotFoundError(f"HiFi FASTQ directory does not exist: {fastq_dir}")
    if config.get("index_bams", False):
        index_bams(config, raw_bams)
    if config.get("convert_bams", False):
        convert_bams(config, raw_bams, fastq_dir)
    if config.get("assemble", True):
        assemble_reads(config, fastq_dir, gfa_dir)
    if config.get("convert_gfa", True):
        if not gfa_dir.is_dir():
            raise FileNotFoundError(f"GFA directory does not exist: {gfa_dir}")
        convert_gfas(config, gfa_dir, fasta_dir)


if __name__ == "__main__":
    main()
