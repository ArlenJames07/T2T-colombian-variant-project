# Structural-variant workflows

Each caller is a standalone Python script with its own private JSON
configuration. The scripts use the same input types and core parameters as the
legacy analyses, write below `results/sv/`, and skip non-empty completed outputs.

## Configure

Create a local configuration for every caller you intend to run:

```bash
cp config/workflows/pbsv.example.json config/workflows/pbsv.local.json
cp config/workflows/cutesv.example.json config/workflows/cutesv.local.json
cp config/workflows/dipcall.example.json config/workflows/dipcall.local.json
cp config/workflows/nanovar.example.json config/workflows/nanovar.local.json
cp config/workflows/svim_asm.example.json config/workflows/svim_asm.local.json
cp config/workflows/sniffles.example.json config/workflows/sniffles.local.json
cp config/workflows/hifi_cnv.example.json config/workflows/hifi_cnv.local.json
cp config/workflows/sawfish.example.json config/workflows/sawfish.local.json
```

Edit each `*.local.json` file with the machine-specific input, reference, and
executable paths. Local configurations are ignored by Git; the example files
are safe to publish.

## Read-based callers

Run each desired caller from the repository root:

```bash
python3 workflow/03_structural_variants/pbsv.py
python3 workflow/03_structural_variants/cuteSV.py
python3 workflow/03_structural_variants/nanovar.py
python3 workflow/03_structural_variants/sniffles.py
```

pbsv performs signature discovery, per-sample calling, and the legacy pass/minimum
length filtering. cuteSV and Sniffles consume aligned BAMs. NanoVar consumes HiFi
FASTQ files.

## Assembly-based callers

Dipcall requires matching haplotype 1 and haplotype 2 FASTAs. SVIM-asm uses the
configured primary assembly suffix:

```bash
python3 workflow/03_structural_variants/dipcall.py
python3 workflow/03_structural_variants/svim_asm.py
```

## Callers that require phased BAMs

Run the phasing workflow first, or point these configurations to an existing
set of phased BAMs and phased small variants:

```bash
python3 workflow/03_structural_variants/hifi_cnv.py
python3 workflow/03_structural_variants/sawfish.py
```

Sawfish runs discovery and joint calling per sample. When
`split_joint_vcfs` is enabled, it also creates separate breakpoint-supported SV
and depth-only CNV VCFs with bcftools and indexes them with tabix.

There are no command-line options. To change a setting, edit only the matching
local configuration, then rerun the same command.
