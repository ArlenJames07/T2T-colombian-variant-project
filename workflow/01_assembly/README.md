# Assembly and scaffolding workflow

This workflow builds phased HiFi assemblies with hifiasm, converts the main GFA
graphs to FASTA, and scaffolds each FASTA against the configured reference with
RagTag. Both scripts resume safely by skipping non-empty final outputs.

## Configuration

Copy and edit the public templates once. The resulting local files are ignored
by Git and may contain absolute paths:

```bash
cp config/workflows/hifiasm.example.json config/workflows/hifiasm.local.json
cp config/workflows/ragtag.example.json config/workflows/ragtag.local.json
```

No command-line arguments are used. Enable or disable the BAM indexing, BAM to
FASTQ conversion, assembly, and GFA conversion stages with the Boolean values in
`hifiasm.local.json`.

## Run the workflow

From the repository root, run assembly first:

```bash
python3 workflow/01_assembly/hifiasm.py
```

The main graphs and FASTAs are written to `results/assemblies/gfa/` and
`results/assemblies/fasta/`.

Then run reference-guided scaffolding:

```bash
python3 workflow/01_assembly/ragtag.py
```

RagTag working outputs are written below
`results/assemblies/scaffolds/<reference>/ragtag/`. Stable relative links to the
scaffold FASTAs are created in the adjacent `fasta/` directory; large FASTA
files are never duplicated by this step.
