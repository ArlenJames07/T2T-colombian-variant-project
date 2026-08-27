# Variant phasing workflow

This workflow uses HiPhase to phase one small-variant VCF and one
structural-variant VCF against each sample's aligned HiFi BAM. It writes phased
BAMs to `results/phasing/bamfiles/` and phased VCFs plus summary tables to
`results/phasing/variants/`.

## Configure

Create and edit the private local configuration:

```bash
cp config/workflows/hiphase.example.json config/workflows/hiphase.local.json
```

Set the directories containing small-variant VCFs, structural-variant VCFs, and
aligned BAMs. Files are matched by the sample code in their names, so the full
input filenames do not have to be identical.

## Run

From the repository root:

```bash
python3 workflow/04_phasing/hiphase.py
```

No command-line arguments are required. A sample is skipped only when its
phased BAM, phased small-variant VCF, and phased structural-variant VCF are all
present and non-empty. Missing input pairs are reported without stopping the
remaining samples.
