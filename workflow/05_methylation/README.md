# Methylation workflow

This workflow runs PacBio `aligned_bam_to_cpg_scores` from pb-CpG-tools on the
phased, haplotagged HiFi BAMs in `results/phasing/bamfiles/`. It produces
combined and haplotype-specific CpG methylation calls in BED and bigWig format
under `results/methylation/<sample>/`.

The default configuration uses the recommended model-based pileup and de novo
CpG-site modes, requires coverage of at least 10 reads, and reads haplotypes
from each alignment's `HP` tag.

## Configure

Create the private local configuration and set the installed executable and
model paths:

```bash
cp config/workflows/pbcpgtools.example.json config/workflows/pbcpgtools.local.json
```

The phased BAM directory already defaults to the output of the phasing stage.
Each BAM must be indexed and must retain the `MM`/`ML` 5mC tags. The `HP` tags
written by HiPhase enable the haplotype-specific outputs.

## Run

From the repository root:

```bash
python3 workflow/05_methylation/pbcpgtools.py
```

The script processes samples sequentially, skips samples whose expected outputs
are already non-empty, and reports phased BAMs that do not yet have a BAI or CSI
index. With `require_haplotype_outputs` enabled, a sample is complete only after
all six combined, haplotype 1, and haplotype 2 BED/bigWig files exist.

Set `modsites_mode` to `reference` and provide `reference_fasta` if calls should
be restricted to CpG sites in the reference rather than the default de novo
sites.
