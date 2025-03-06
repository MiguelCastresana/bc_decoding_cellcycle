# --- Install helper packages if needed ---
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# --- Set a CRAN mirror (if not already set) ---
chooseCRANmirror(graphics = FALSE)

# --- Install specific versions of CRAN packages ---
# Note: For some packages the CRAN name may differ from the conda name.
# For example, r-seurat is available as "Seurat" on CRAN.
remotes::install_version("tidyverse", version = "2.0.0", upgrade = "never")
remotes::install_version("Seurat", version = "5.2.0", upgrade = "never")
remotes::install_version("SeuratObject", version = "5.0.2", upgrade = "never")
remotes::install_version("Matrix", version = "1.6-4", upgrade = "never")      # r-matrix
remotes::install_version("readxl", version = "1.4.3", upgrade = "never")
remotes::install_version("yaml", version = "2.3.10", upgrade = "never")
remotes::install_version("msigdbr", version = "7.5.1", upgrade = "never")
remotes::install_version("pheatmap", version = "1.0.12", upgrade = "never")
remotes::install_version("caret", version = "7.0-1", upgrade = "never")
remotes::install_version("devtools", version = "2.4.5", upgrade = "never")
remotes::install_version("igraph", version = "2.1.3", upgrade = "never")
# dplyr, forcats, ggplot2, readr, tidyr, tibble are included in tidyverse, but you can force reinstall if needed:
remotes::install_version("dplyr", version = "1.1.4", upgrade = "never")
remotes::install_version("forcats", version = "1.0.0", upgrade = "never")
remotes::install_version("ggplot2", version = "3.5.1", upgrade = "never")
remotes::install_version("readr", version = "2.1.5", upgrade = "never")
remotes::install_version("tidyr", version = "1.3.1", upgrade = "never")
remotes::install_version("tibble", version = "3.2.1", upgrade = "never")
install.packages("scCustomize", type = "binary")

# --- Install GitHub packages at specific tags ---
remotes::install_github("aertslab/SCopeLoomR",force = TRUE)
remotes::install_github("stephenturner/annotables", ref = "v0.2.0", force = TRUE)
remotes::install_github("yanlinlin82/ggvenn",version = "0.1.16", force = TRUE)
remotes::install_github("aertslab/SCENIC", version = "1.3.1", force = TRUE)



# --- Install Bioconductor packages ---
# Set the Bioconductor release that corresponds to your desired versions.
# For example, if you expect:
#   singler ~2.8.0, SummarizedExperiment ~1.32.0, scuttle ~1.16.0,
#   edgeR ~4.4.1, clusterProfiler ~4.14.4, SingleCellExperiment ~1.28.1, 
#   and org.Hs.eg.db ~3.20.0,
# then choose a Bioconductor release (e.g., "3.2") known to provide those.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
BiocManager::install("singler")                    # bioconductor-singler (expected ~2.8.0)
BiocManager::install("SummarizedExperiment")         # bioconductor-summarizedexperiment (expected ~1.32.0)
BiocManager::install("scuttle")                      # bioconductor-scuttle (expected ~1.16.0)
BiocManager::install("edgeR")                        # bioconductor-edger (expected ~4.4.1)
BiocManager::install("clusterProfiler")              # bioconductor-clusterprofiler (expected ~4.14.4)
BiocManager::install("SingleCellExperiment")         # bioconductor-singlecellexperiment (expected ~1.28.1)
BiocManager::install("org.Hs.eg.db")                 # bioconductor-org.hs.eg.db (expected ~3.20.0)
BiocManager::install("limma")                
BiocManager::install("fgsea")  
