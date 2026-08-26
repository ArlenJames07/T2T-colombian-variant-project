# T2T Colombian Variant Project

Reproducible code for PacBio HiFi assembly, small-variant and
structural-variant calling, phasing, methylation, annotation, and multi-omic
integration against T2T-CHM13v2.0.

The expensive analyses have already been run. The directories under `results/`
are symbolic links to those existing outputs; no data were copied and no caller
was rerun during the reorganization. See [docs/pipeline.md](docs/pipeline.md) for
the exact link map and workflow boundaries.

## Layout

```text
config/       cohort inputs, paths, and parameters
metadata/     phenotype and file manifests
workflow/     one import-safe CLI module per pipeline stage
analyses/     paper-oriented analysis workspaces
resources/    references and external annotations
results/      links to completed primary and integrated results
figures/      links to existing figure collections
reports/      links to finalized reports and tables
docs/         pipeline, data dictionary, and analysis plan
```

## Running a module

Every workflow file is a standalone CLI and does nothing when imported. Use
`--help` before running it. Most compute-heavy modules support `--dry-run`,
`--force`, and explicit `--output` paths.

```bash
python3 workflow/03_structural_variants/sniffles.py --help
python3 workflow/04_phasing/hiphase.py --help
python3 workflow/05_methylation/pbcpgtools.py --help
```

## HiFi assembly workflow

The assembly workflow runs hifiasm and converts its primary-contig GFA outputs
to FASTA. It is resumable: existing, non-empty outputs are skipped by default,
so completed samples are not assembled again.

The completed assembly data remain on the large `/mnt` disk and are exposed as
short, per-file symbolic links under stable repository directories:

```text
results/assemblies/
├── fasta/
│   ├── 001P_hap1.fasta
│   ├── 001P_hap2.fasta
│   └── 001P_p_ctg.fasta
└── gfa/
    ├── 001P_hap1.gfa
    ├── 001P_hap2.gfa
    └── 001P_p_ctg.gfa
```

Each name contains only the cohort sample code and assembly type. The linked
FASTA sources are in `/mnt/diskrare/arlenb/03_fasta_files`, and the linked GFA
sources are in `/mnt/diskrare/arlenb/02_hifi_results`.

From the repository root, inspect the planned work first:

```bash
cd /home/rare/arlen/T2T-colombian-variant-project
python3 workflow/01_assembly/hifiasm.py --dry-run
```

For the current completed dataset, the summary should report zero commands to
run and 72 existing outputs skipped. To resume and create only missing outputs:

```bash
python3 workflow/01_assembly/hifiasm.py
```

Individual stages or samples can be selected explicitly:

```bash
# Run or inspect only hifiasm
python3 workflow/01_assembly/hifiasm.py --step assemble --dry-run

# Run or inspect only GFA-to-FASTA conversion
python3 workflow/01_assembly/hifiasm.py --step convert --dry-run

# Restrict processing to a sample-name substring
python3 workflow/01_assembly/hifiasm.py --sample 001P --dry-run

# Show every available option
python3 workflow/01_assembly/hifiasm.py --help
```

Do not use `--force` unless completed outputs should intentionally be replaced.
The withdrawn alternate run `bc2044v2` is listed under
`assembly.exclude_name_patterns` in `config/params.yaml`, which prevents resume
mode from recreating its removed FASTA files.

The repository paths are used by the workflow and documentation, while the
large files remain on `/mnt`; no data need to be copied into Git. On another
machine, update `config/paths.yaml` and recreate the result links for the local
storage layout.

Machine-specific locations are centralized in `config/paths.yaml`; algorithmic
settings are in `config/params.yaml`. The original source scripts remain in
`/home/rare/arlen/scripts` and `/home/rare/arlen/genome_assembly.py` as a
provenance record.
