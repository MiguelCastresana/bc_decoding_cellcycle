# Data Notes

This repository stores analysis code, configuration, and documentation. Raw single-cell atlases, intermediate objects, and generated figures are intentionally kept outside Git because they are large and often have separate access requirements.

## External Inputs

Expected input categories:

- Discovery breast cancer single-cell atlas.
- Validation breast cancer single-cell atlas.
- Cell type annotation references for SingleR.
- PAM50/subtype marker resources.
- Pathway resources used by enrichment scripts.
- SCENIC databases and motif resources used by the GRN notebook.
- DrugBank and CTD resources for drug-prioritization steps.

## Local Layout

Use a local data directory that is not committed, for example:

```text
data/
├── raw/
├── processed/
├── scenic/
└── drug_resources/
```

If a script expects a different path, update the corresponding YAML config under `src/DEG/` or `src/GRN/`.

## Git Hygiene

Generated tables, R objects, figures, and local session files are ignored by `.gitignore`. Commit small example inputs only when they are explicitly meant for tests or documentation.
