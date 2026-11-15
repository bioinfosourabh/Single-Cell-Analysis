#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R -- run the numbered scripts in order.
#
# Usage (from the repository root):
#   Rscript run_pipeline.R                                   # PBMC 3k (config.yaml)
#   Rscript run_pipeline.R --config config_breast_tumour.yaml
#   Rscript run_pipeline.R --from 04 --to 06
#
# Step 00 downloads the data if it is not already present. Each step runs in a
# fresh R process; the run stops at the first failure.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) > i) return(args[i + 1])
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", flag, "="), "", hit[1]))
  default
}

config <- get_arg("--config", "config.yaml")
from <- get_arg("--from", NULL)
to <- get_arg("--to", NULL)
if (!file.exists(config)) stop("Configuration not found: ", config)

steps <- c("00_download_data.R", "01_qc.R", "02_normalize_pca.R", "03_cluster_umap.R",
           "04_markers.R", "05_annotation_singler.R", "06_subcluster_analysis.R")
ids <- sub("_.*$", "", steps)
if (!is.null(from)) steps <- steps[seq(match(from, ids), length(steps))]
ids <- sub("_.*$", "", steps)
if (!is.null(to)) steps <- steps[seq_len(match(to, ids))]

rscript <- file.path(R.home("bin"), "Rscript")
t0 <- Sys.time()
for (s in steps) {
  cat("\n>>> ", s, "\n", sep = "")
  t1 <- Sys.time()
  status <- system2(rscript, c(file.path("scripts", s), "--config", config))
  if (status != 0) stop("Step failed: ", s, " (exit status ", status, ")")
  cat(sprintf(">>> %s finished in %.1f s\n", s, as.numeric(difftime(Sys.time(), t1, units = "secs"))))
}
cat(sprintf("\nPipeline finished in %.1f min (%d steps, config: %s)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), length(steps), config))
