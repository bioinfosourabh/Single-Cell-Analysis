#!/usr/bin/env Rscript
# =============================================================================
# 02_normalize_pca.R
#
# Normalisation (LogNormalize, or SCTransform v2 if normalization$method is
# "SCT"), highly variable genes, scaling and PCA.
#
# Output: results/rds/02_seurat_pca.rds
#         results/tables/top_variable_features.csv
#         results/figures/Highly_Variable_Features.png, PCA_loadings.png,
#         PCA_Plot.png, Dimheatmap_for_first_2_PC.png, ElbowPlot_PCA.png
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

cfg <- load_config()
nc <- cfg$normalization
log_step(paste0("02. Normalisation (", nc$method, "), variable features and PCA"))

require_files(p_rds(cfg, "01_seurat_qc.rds"))
obj <- readRDS(p_rds(cfg, "01_seurat_qc.rds"))

if (identical(nc$method, "SCT")) {
  obj <- SCTransform(obj, vst.flavor = "v2", variable.features.n = nc$n_variable_features,
                     vars.to.regress = "percent.mt", verbose = FALSE)
  # Log-normalised RNA data are still needed for marker plots and SingleR.
  obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                       scale.factor = nc$scale_factor, verbose = FALSE)
  hvg_plot <- NULL
} else {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = nc$scale_factor, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = nc$n_variable_features, verbose = FALSE)
  # Scale all genes so that heatmaps of arbitrary marker genes are possible.
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
  hvg_plot <- TRUE
}

hvg <- VariableFeatures(obj)
top10 <- head(hvg, 10)
log_msg("Variable features: ", length(hvg), "; top 10: ", paste(top10, collapse = ", "))
utils::write.csv(data.frame(rank = seq_along(hvg), gene = hvg),
                 p_table(cfg, "top_variable_features.csv"), row.names = FALSE)

if (!is.null(hvg_plot)) {
  p <- LabelPoints(plot = VariableFeaturePlot(obj), points = top10, repel = TRUE, xnudge = 0, ynudge = 0) +
    ggtitle("Highly variable features (vst)") + theme_minimal()
  save_plot(cfg, p, "Highly_Variable_Features", 9, 6)
}

n_pcs <- min(50, ncol(obj) - 1, length(hvg) - 1)
obj <- RunPCA(obj, features = hvg, npcs = n_pcs, verbose = FALSE)
print(obj[["pca"]], dims = 1:min(5, n_pcs), nfeatures = 5)

save_plot(cfg, VizDimLoadings(obj, dims = 1:min(6, n_pcs), nfeatures = 15, reduction = "pca",
                              ncol = 3),
          "PCA_loadings", 12, 10)
save_plot(cfg, DimPlot(obj, reduction = "pca") + ggtitle("PCA") + NoLegend(),
          "PCA_Plot", 7, 6)
save_base_plot(cfg, function() {
  DimHeatmap(obj, dims = 1:2, cells = min(500, ncol(obj)), balanced = TRUE, fast = TRUE)
}, "Dimheatmap_for_first_2_PC", 10, 6)
save_plot(cfg, ElbowPlot(obj, ndims = n_pcs) + ggtitle("Elbow plot"), "ElbowPlot_PCA", 7, 5)

saveRDS(obj, p_rds(cfg, "02_seurat_pca.rds"))
write_session_info(cfg, "02_normalize_pca")
