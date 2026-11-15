#!/usr/bin/env Rscript
# =============================================================================
# 00_download_data.R
#
# Download and unpack the 10x Genomics PBMC 3k filtered gene-barcode matrices
# (the dataset of the Seurat guided clustering tutorial) into
# data/pbmc3k/filtered_gene_bc_matrices/hg19/.
#
# Usage: Rscript scripts/00_download_data.R [--config config.yaml]
# =============================================================================

source("scripts/utils.R")

cfg <- load_config()
log_step("00. Download data")

mtx <- file.path(cfg$paths$data_dir, c("matrix.mtx", "matrix.mtx.gz"))
if (any(file.exists(mtx))) {
  log_msg("Data already present in ", cfg$paths$data_dir, "; nothing to do.")
} else if (is.null(cfg$download$url)) {
  stop("No data in ", cfg$paths$data_dir, " and no download URL in ",
       basename(cfg$config_file), ". Place the 10x matrix there manually.")
} else {
  options(timeout = max(600, getOption("timeout")))
  dir.create(cfg$download$dest_dir, recursive = TRUE, showWarnings = FALSE)
  archive <- file.path(cfg$download$dest_dir, basename(cfg$download$url))
  log_msg("Downloading ", cfg$download$url)
  utils::download.file(cfg$download$url, archive, mode = "wb", quiet = TRUE)
  if (!is.null(cfg$download$md5)) {
    md5 <- unname(tools::md5sum(archive))
    if (!identical(md5, cfg$download$md5)) stop("MD5 mismatch for ", archive, ": ", md5)
  }
  utils::untar(archive, exdir = cfg$download$dest_dir)
  unlink(archive)
  if (!any(file.exists(mtx))) {
    stop("matrix.mtx not found in ", cfg$paths$data_dir, " after extraction; ",
         "check paths$data_dir in the configuration.")
  }
  log_msg("Extracted to ", cfg$paths$data_dir)
}
print(list.files(cfg$paths$data_dir))

write_session_info(cfg, "00_download_data")
