#!/usr/bin/env python3
"""Call assembly-based variants from paired haplotype assemblies with dipcall."""

import json
import re
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "workflows" / "dipcall.local.json"


def load_config():
    if not CONFIG_FILE.is_file():
        raise FileNotFoundError(f"Missing {CONFIG_FILE}; copy dipcall.example.json to dipcall.local.json.")
    with CONFIG_FILE.open(encoding="utf-8") as handle:
        return json.load(handle)


def path_value(config, key):
    value = Path(config[key]).expanduser()
    return value if value.is_absolute() else PROJECT_ROOT / value


def sample_id(name):
    match = re.search(r"(\d{3}[A-Za-z])", name)
    return match.group(1).upper() if match else name.split("_")[0]


def main():
    config = load_config()
    assemblies = path_value(config, "assemblies_dir")
    reference = path_value(config, "reference_fasta")
    output_dir = path_value(config, "output_dir") / config.get("reference_name", "t2t")
    if not assemblies.is_dir():
        raise FileNotFoundError(f"Assembly directory does not exist: {assemblies}")
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
    output_dir.mkdir(parents=True, exist_ok=True)
    hap1_suffix = config.get("haplotype_1_suffix", "hap1.fa")
    hap2_suffix = config.get("haplotype_2_suffix", "hap2.fa")

    for hap1 in sorted(assemblies.glob(f"*{hap1_suffix}")):
        hap2 = Path(str(hap1).replace(hap1_suffix, hap2_suffix))
        if not hap2.is_file():
            raise FileNotFoundError(f"Missing haplotype 2 assembly for {hap1.name}: {hap2}")
        sample = sample_id(hap1.name)
        prefix = output_dir / sample
        final_vcf = output_dir / f"{sample}.dip.vcf"
        if final_vcf.is_file() and final_vcf.stat().st_size > 0:
            print(f"Skipping completed sample: {sample}")
            continue
        makefile = output_dir / f"{sample}.mak"
        temporary = makefile.with_suffix(".mak.tmp")
        command = [config["dipcall"], prefix, reference, hap1, hap2, "-t", config.get("threads", 32)]
        print("Running:", " ".join(map(str, command)), flush=True)
        with temporary.open("wb") as output:
            subprocess.run([str(value) for value in command], check=True, stdout=output)
        temporary.replace(makefile)
        subprocess.run(["make", f"-j{config.get('make_jobs', 2)}", "-f", str(makefile)], check=True)


if __name__ == "__main__":
    main()
