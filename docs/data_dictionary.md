# Data dictionary

## `config/samples.tsv`

| Column | Meaning |
|---|---|
| `sample_id` | Stable short cohort identifier |
| `bam` | Original PacBio BAM filename |
| `group` | PWS, AS, DiGeorge, or Control |
| `mechanism` | Molecular mechanism or Control |
| `sex` | Reported sex (`F` or `M`) |
| `age` | Age at sampling, in years |
| `include` | Whether workflow batch discovery should include the sample |

## Integration tables

The generic integration modules use these canonical fields:

| Field | Meaning |
|---|---|
| `sample_id` | Sample identifier |
| `sv_id` | Stable structural-variant identifier |
| `chrom`, `start`, `end` | T2T genomic interval |
| `breakpoint`, `side` | Breakpoint coordinate and left/right designation |
| `haplotype` | Haplotype label or `unassigned` |
| `cpg_sites` | CpGs contributing to a summary |
| `mean_methylation` | Mean methylation score in the interval |
| `median_methylation` | Median methylation score in the interval |

Coordinates should be documented as zero- or one-based by every producing
analysis. `breakpoint_windows.py` preserves input coordinates and clips negative
window starts to zero.
