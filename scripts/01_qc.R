#!/usr/bin/env Rscript
# =============================================================================
# 01_qc.R
#
# Load the 10x matrix, create the Seurat object, compute mitochondrial and
# ribosomal read fractions, plot QC metrics and filter cells.
#
# Output: results/rds/01_seurat_qc.rds
#         results/tables/qc_cell_metrics_prefilter.csv, qc_filter_summary.csv
#         results/figures/Pre_filter.png, QC_Metric_relationships.png,
#         After_filtering.png
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

cfg <- load_config()
q <- cfg$qc
log_step("01. Load data and quality control")

counts <- Read10X(data.dir = cfg$paths$data_dir)
if (is.list(counts)) counts <- counts[["Gene Expression"]]
rownames(counts) <- gsub("_", "-", rownames(counts))
log_msg("Raw matrix: ", nrow(counts), " genes x ", ncol(counts), " barcodes")

obj <- CreateSeuratObject(counts = counts, project = cfg$project,
                          min.cells = q$min_cells, min.features = q$min_features)
rm(counts)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = q$mt_pattern)
obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = q$rb_pattern)
log_msg("After min.cells / min.features: ", nrow(obj), " genes x ", ncol(obj), " cells")
if (sum(grepl(q$mt_pattern, rownames(obj))) == 0) {
  warning("No genes match the mitochondrial pattern '", q$mt_pattern, "'.")
}

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
utils::write.csv(obj@meta.data, p_table(cfg, "qc_cell_metrics_prefilter.csv"))

vln <- function(o, title) {
  VlnPlot(o, features = qc_features, ncol = 4, pt.size = 0.1, layer = "counts") +
    plot_annotation(title = title)
}
save_plot(cfg, vln(obj, "QC metrics before filtering"), "Pre_filter", 16, 4.5)

p1 <- FeatureScatter(obj, "nCount_RNA", "percent.mt") + ggtitle("UMI counts vs % mitochondrial")
p2 <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA") + ggtitle("UMI counts vs genes detected")
p3 <- FeatureScatter(obj, "nCount_RNA", "percent.rb") + ggtitle("UMI counts vs % ribosomal")
save_plot(cfg, (p1 | p2 | p3) & theme(legend.position = "none"),
          "QC_Metric_relationships", 16, 5)

# ---- Filtering -----------------------------------------------------------------
md <- obj@meta.data
keep <- md$nFeature_RNA >= q$min_features & md$percent.mt < q$max_percent_mt
if (!is.null(q$max_features)) keep <- keep & md$nFeature_RNA < q$max_features
if (!is.null(q$max_counts)) keep <- keep & md$nCount_RNA < q$max_counts
if (!is.null(q$max_percent_rb)) keep <- keep & md$percent.rb < q$max_percent_rb

n_before <- ncol(obj)
obj <- subset(obj, cells = colnames(obj)[keep])
summary_df <- data.frame(
  step = c("cells before filtering", "cells after filtering", "cells removed"),
  n = c(n_before, ncol(obj), n_before - ncol(obj)))
print(summary_df, row.names = FALSE)
utils::write.csv(summary_df, p_table(cfg, "qc_filter_summary.csv"), row.names = FALSE)

save_plot(cfg, vln(obj, "QC metrics after filtering"), "After_filtering", 16, 4.5)
saveRDS(obj, p_rds(cfg, "01_seurat_qc.rds"))

write_session_info(cfg, "01_qc")
