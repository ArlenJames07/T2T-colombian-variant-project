# Workflow configuration

Each workflow reads its matching `*.local.json` file automatically. Local files
contain machine-specific absolute paths and are ignored by Git. The
`*.example.json` files are safe templates intended for the remote repository.

To configure a workflow on another machine, copy its example file to the same
name with `.example.json` replaced by `.local.json`, then edit the paths. No
command-line arguments are required or accepted:

```bash
python3 workflow/01_assembly/hifiasm.py
```

All relative paths are resolved from the repository root. Every script skips a
sample when its final non-empty output already exists.
