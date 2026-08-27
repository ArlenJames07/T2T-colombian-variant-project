# Small-variant workflow

This workflow calls SNVs and small indels from reference-aligned PacBio HiFi
reads with DeepVariant. One compressed VCF and one compressed gVCF are produced
per sample under `results/snvs/<reference>/`.

## Configure

Create the private local configuration from the public template and edit the
input, reference, Docker image, and thread settings:

```bash
cp config/workflows/deepvariant.example.json config/workflows/deepvariant.local.json
```

The BAM directory and reference must be accessible to the Docker daemon. The
script mounts both inputs read-only and mounts only the project output directory
as writable.

## Run

From the repository root:

```bash
python3 workflow/02_small_variants/deepvariant.py
```

No arguments are required. A sample with a non-empty final VCF already present
is skipped, so the same command can resume an interrupted run without replacing
completed calls.
