#!/usr/bin/env Rscript
# =============================================================================
# install_packages.R -- install the R packages used by the workflow.
# Only missing packages are installed. Usage: Rscript install_packages.R
# =============================================================================

cran_repo <- unname(getOption("repos")["CRAN"])
if (length(cran_repo) == 0 || is.na(cran_repo) || cran_repo == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
ncpus <- max(1L, parallel::detectCores() - 1L)

cran <- c("Seurat", "SeuratObject", "Matrix", "dplyr", "ggplot2", "patchwork",
          "RColorBrewer", "yaml", "remotes", "BiocManager")
bioc <- c("SingleR", "celldex", "SingleCellExperiment", "SummarizedExperiment",
          "glmGamPoi")   # glmGamPoi speeds up SCTransform v2

missing <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install <- missing(c(cran, bioc))
if (length(to_install) > 0) {
  message("Installing: ", paste(to_install, collapse = ", "))
  BiocManager::install(to_install, update = FALSE, ask = FALSE, Ncpus = ncpus)
}

# Optional: presto provides a fast Wilcoxon test that Seurat uses automatically.
if (!requireNamespace("presto", quietly = TRUE)) {
  tryCatch(remotes::install_github("immunogenomics/presto", upgrade = "never"),
           error = function(e) message("presto not installed (optional): ", conditionMessage(e)))
}

still_missing <- missing(c(cran, bioc))
if (length(still_missing) > 0) {
  stop("Packages failed to install: ", paste(still_missing, collapse = ", "))
}
message("All packages available. R ", getRversion(), ", Seurat ", packageVersion("Seurat"))
