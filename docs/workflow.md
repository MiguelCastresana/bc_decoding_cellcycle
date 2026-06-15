# Workflow

Run scripts from the repository root unless a script or notebook states otherwise.

## 1. Data Preparation

- `src/data/data_generation_discovery.R`
- `src/data/data_generation_validation.R`
- `src/data/common_genes.R`

These scripts prepare discovery and validation single-cell objects, normalize data, select features, and align the gene universe used downstream.

## 2. Cell Annotation And Subtyping

- `src/data/singleR_training_celltype.R`
- `src/data/singleR_validation.R`
- `src/data/sc_subtyping_discovery.R`
- `src/data/sc_subtyping_validation.R`

This stage assigns cell types and breast cancer subtype labels for downstream phase-subtype comparisons.

## 3. Differential Expression And Pathways

- `src/DEG/DGE_analysis.R`
- `src/DEG/intersection.R`
- `src/DEG/DE_across_figure3.R`
- `src/DEG/DE_within_figure3.R`
- `src/DEG/DE_venn_diagram.R`

Configuration files in `src/DEG/` define input/output locations for these analyses.

## 4. Gene Regulatory Networks

- `src/GRN/GRN_analysis.ipynb`
- `src/GRN/regulon_analysis.R`
- `src/GRN/regulon_analysis_nophase.R`
- `src/GRN/intersection.R`

The Python notebook handles SCENIC-style GRN inference. The R scripts summarize regulons, pathway signals, and intersections across datasets.

## 5. Figures

Figure scripts live near the analysis stage they visualize:

- `src/DEG/DE_*figure3.R`
- `src/GRN/GRN_figure4*.R`
- `src/GRN/GRN_figure5.py`

Keep generated figure files outside Git unless they are intentionally used in the README.
