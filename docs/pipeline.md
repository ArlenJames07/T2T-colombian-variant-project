# Pipeline organization

The rewritten scripts in `workflow/` do not accept command-line arguments.
Each script reads a private `config/workflows/*.local.json` file containing the
input and executable paths for the current machine. Git ignores local files;
only their `*.example.json` templates should be published.

| Stage | Script | Local configuration | Project output |
|---|---|---|---|
| Assembly | `workflow/01_assembly/hifiasm.py` | `hifiasm.local.json` | `results/assemblies/` |
| Scaffolding | `workflow/01_assembly/ragtag.py` | `ragtag.local.json` | `results/assemblies/scaffolds/` |
| Small variants | `workflow/02_small_variants/deepvariant.py` | `deepvariant.local.json` | `results/snvs/` |
| pbsv | `workflow/03_structural_variants/pbsv.py` | `pbsv.local.json` | `results/sv/pbsv/` |
| HiFiCNV | `workflow/03_structural_variants/hifi_cnv.py` | `hifi_cnv.local.json` | `results/sv/hifi_cnv/` |
| cuteSV | `workflow/03_structural_variants/cuteSV.py` | `cutesv.local.json` | `results/sv/cutesv/` |
| dipcall | `workflow/03_structural_variants/dipcall.py` | `dipcall.local.json` | `results/sv/dipcall/` |
| NanoVar | `workflow/03_structural_variants/nanovar.py` | `nanovar.local.json` | `results/sv/nanovar/` |
| Sawfish | `workflow/03_structural_variants/sawfish.py` | `sawfish.local.json` | `results/sv/sawfish/` |
| SVIM-asm | `workflow/03_structural_variants/svim_asm.py` | `svim_asm.local.json` | `results/sv/svim_asm/` |
| Sniffles | `workflow/03_structural_variants/sniffles.py` | `sniffles.local.json` | `results/sv/sniffles/` |
| Phasing | `workflow/04_phasing/hiphase.py` | `hiphase.local.json` | `results/phasing/` |
| Methylation | `workflow/05_methylation/pbcpgtools.py` | `pbcpgtools.local.json` | `results/methylation/` |

Run any configured workflow directly from the repository root, for example:

```bash
python3 workflow/01_assembly/hifiasm.py
```

Every script creates its output directory and skips a sample when the expected
final output already exists and is non-empty. This behavior makes reruns safe
for completed or partially completed datasets.

Previously generated results remain in their original storage locations and
are exposed locally through links where needed. They are not recopied during
repository organization. The `results/` tree and all private local
configurations are ignored by Git, so neither large data nor local absolute
paths are included in the remote repository.
