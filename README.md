# T2T Colombian Variant Project

Reproducible workflows for PacBio HiFi assembly, small-variant calling,
structural-variant calling, and phasing against T2T-CHM13v2.0.

The compute-intensive analyses were completed before this reorganization.
Existing results are referenced locally under `results/`; this repository does
not copy those large files and Git ignores the complete results tree. The new
scripts can independently recreate missing outputs in the same project folders.

## Repository layout

```text
config/workflows/         public examples and private per-workflow settings
workflow/01_assembly/     hifiasm and RagTag
workflow/02_small_variants/
workflow/03_structural_variants/
workflow/04_phasing/
results/                  generated or locally linked results (Git-ignored)
metadata/                 sample and cohort tables
resources/                references and external resources
docs/                     pipeline and analysis documentation
```

## Configuration

Every script automatically reads one matching `config/workflows/*.local.json`
file. These local files contain the original machine-specific input paths and
are ignored by Git. The corresponding `*.example.json` files contain only safe
placeholder paths and are intended for the remote repository.

For example:

```bash
cp config/workflows/hifiasm.example.json config/workflows/hifiasm.local.json
# Edit hifiasm.local.json, then run:
python3 workflow/01_assembly/hifiasm.py
```

No workflow accepts or requires command-line arguments. Relative paths in the
JSON configuration are resolved from the repository root. Each workflow skips
non-empty final outputs, making interrupted runs resumable without overwriting
completed samples.

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
```

Do not add large outputs or local result links to Git. Because all of
`results/` except its `.gitkeep` file is ignored, link targets and private local
storage paths are not exposed in the remote repository.
