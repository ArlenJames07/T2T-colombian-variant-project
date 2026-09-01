#!/usr/bin/env python3
"""Call CpG methylation from phased, haplotagged HiFi BAMs."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "pbcpgtools.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(
            f"Missing {CONFIG_FILE}; copy pbcpgtools.example.json to "
            "pbcpgtools.local.json."
        )
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    value = Path(config[key]).expanduser()
    return value if value.is_absolute() else PROJECT_ROOT / value


def sample_id(name):
    match = re.search(r"(\d{3}[A-Za-z])(?:\D|$)", name)
    return match.group(1).upper() if match else name.split(".")[0]


def complete(path):
    return path.is_file() and path.stat().st_size > 0


def bam_index(bam):
    """Return an existing BAI or CSI index for a BAM, if one is available."""
    candidates = (
        Path(f"{bam}.bai"),
        bam.with_suffix(".bai"),
        Path(f"{bam}.csi"),
        bam.with_suffix(".csi"),
    )
    return next((path for path in candidates if complete(path)), None)


def expected_outputs(prefix, require_haplotype_outputs):
    outputs = [
        Path(f"{prefix}.combined.bed"),
        Path(f"{prefix}.combined.bw"),
    ]
    if require_haplotype_outputs:
        for haplotype in ("hap1", "hap2"):
            outputs.extend(
                [Path(f"{prefix}.{haplotype}.bed"), Path(f"{prefix}.{haplotype}.bw")]
            )
    return outputs


def validate_config(config):
    pileup_mode = config.get("pileup_mode", "model")
    modsites_mode = config.get("modsites_mode", "denovo")
    if pileup_mode not in {"model", "count"}:
        raise ValueError("pileup_mode must be 'model' or 'count'")
    if modsites_mode not in {"denovo", "reference"}:
        raise ValueError("modsites_mode must be 'denovo' or 'reference'")
    defaults = {"threads": 32, "minimum_coverage": 10, "minimum_mapq": 1}
    for key, default in defaults.items():
        if int(config.get(key, default)) < 1:
            raise ValueError(f"{key} must be at least 1")


def main():
    config = load_config()
    validate_config(config)

    bams = path_value(config, "phased_bams")
    output_root = path_value(config, "output_dir")
    executable = path_value(config, "pbcpgtools")
    pileup_mode = config.get("pileup_mode", "model")
    modsites_mode = config.get("modsites_mode", "denovo")

    if not bams.is_dir():
        raise FileNotFoundError(f"Phased BAM directory does not exist: {bams}")
    if not executable.is_file():
        raise FileNotFoundError(f"pb-CpG-tools executable does not exist: {executable}")

    model = None
    if pileup_mode == "model":
        model = path_value(config, "model")
        if not model.is_file():
            raise FileNotFoundError(f"pb-CpG-tools model does not exist: {model}")

    reference = None
    if modsites_mode == "reference":
        reference = path_value(config, "reference_fasta")
        if not reference.is_file():
            raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")

    output_root.mkdir(parents=True, exist_ok=True)
    bam_suffix = config.get("bam_suffix", ".bam")
    input_bams = sorted(bams.glob(f"*{bam_suffix}"))
    if not input_bams:
        raise FileNotFoundError(f"No *{bam_suffix} files found in {bams}")

    require_haplotypes = config.get("require_haplotype_outputs", True)
    for bam in input_bams:
        sample = sample_id(bam.name)
        if bam_index(bam) is None:
            print(f"Skipping {sample}: BAM index not found for {bam}")
            continue

        sample_dir = output_root / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        prefix = sample_dir / sample
        outputs = expected_outputs(prefix, require_haplotypes)
        if all(complete(path) for path in outputs):
            print(f"Skipping completed sample: {sample}")
            continue

        command = [
            executable,
            "--bam", bam,
            "--output-prefix", prefix,
            "--pileup-mode", pileup_mode,
            "--modsites-mode", modsites_mode,
            "--threads", config.get("threads", 32),
            "--min-coverage", config.get("minimum_coverage", 10),
            "--min-mapq", config.get("minimum_mapq", 1),
            "--hap-tag", config.get("haplotype_tag", "HP"),
        ]
        if model is not None:
            command.extend(["--model", model])
        if reference is not None:
            command.extend(["--ref", reference])

        print("Running:", " ".join(map(str, command)), flush=True)
        subprocess.run([str(value) for value in command], check=True)

        missing_outputs = [path for path in outputs if not complete(path)]
        if missing_outputs:
            missing = ", ".join(str(path) for path in missing_outputs)
            raise RuntimeError(f"pb-CpG-tools did not create expected outputs: {missing}")


if __name__ == "__main__":
    main()
