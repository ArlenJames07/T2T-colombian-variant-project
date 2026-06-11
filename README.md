# Genome Assembly and Variant Discovery Pipeline

This repository contains a reproducible workflow for PacBio HiFi genome
assembly, read-based variant calling, assembly-based variant discovery,
variant harmonization, consensus merging, annotation, and downstream analysis.

## Workflow

```text
PacBio HiFi reads
        │
        ├── Read quality control
        │
        ├── Genome assembly
        │      ├── Primary assembly
        │      ├── Haplotype 1
        │      └── Haplotype 2
        │
        ├── Read-to-reference alignment
        │
        ├── Read-based variant calling
        │      ├── Small variants
        │      └── Structural variants
        │
        ├── Assembly-to-reference alignment
        │
        ├── Assembly-based variant calling
        │
        ├── Variant harmonization
        │
        ├── Callset merging
        │      ├── Union callset
        │      └── High-confidence consensus callset
        │
        ├── Variant annotation
        │
        └── Downstream analysis