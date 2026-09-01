# T2T Colombian Variant Project

Reproducible workflows for PacBio HiFi assembly, small-variant calling,
structural-variant calling, phasing, and CpG methylation calling against
T2T-CHM13v2.0.


## Repository layout

```text
config/workflows/         public examples and private per-workflow settings
workflow/01_assembly/     hifiasm and RagTag
workflow/02_small_variants/
workflow/03_structural_variants/
workflow/04_phasing/
workflow/05_methylation/
results/                  generated results 
metadata/                 sample and cohort tables
resources/                references and external resources
docs/                     pipeline and analysis documentation
```

## Configuration

Every script automatically reads one matching `config/workflows/*.local.json`
file. 

For example:

```bash
cp config/workflows/hifiasm.example.json config/workflows/hifiasm.local.json
# Edit hifiasm.local.json, then run:
python3 workflow/01_assembly/hifiasm.py
```

## Run the workflows

Run commands from the repository root. Follow each directory's README for
required inputs and detailed ordering.

Assembly and scaffolding:

```bash
python3 workflow/01_assembly/hifiasm.py
python3 workflow/01_assembly/ragtag.py
```

Small variants:

```bash
python3 workflow/02_small_variants/deepvariant.py
```

Structural variants:

```bash
python3 workflow/03_structural_variants/pbsv.py
python3 workflow/03_structural_variants/cuteSV.py
python3 workflow/03_structural_variants/dipcall.py
python3 workflow/03_structural_variants/nanovar.py
python3 workflow/03_structural_variants/svim_asm.py
python3 workflow/03_structural_variants/sniffles.py
python3 workflow/03_structural_variants/hifi_cnv.py
python3 workflow/03_structural_variants/sawfish.py
```

Phasing:

```bash
python3 workflow/04_phasing/hiphase.py
```

Methylation calling from phased BAMs:

```bash
python3 workflow/05_methylation/pbcpgtools.py
```

The structural-variant README distinguishes callers that use aligned reads,
assemblies, or phased BAMs. HiFiCNV and Sawfish should be run only after the
required phased inputs exist.

## Results and privacy

All new output paths stay below:

```text
results/assemblies/
results/snvs/
results/sv/
results/phasing/
results/methylation/
```

