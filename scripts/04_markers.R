#!/usr/bin/env Rscript
# =============================================================================
# 04_markers.R
#
# Cluster marker genes (Wilcoxon rank-sum test on log-normalised RNA data,
# Bonferroni-adjusted p-values as reported by Seurat) and marker plots.
#
# Output: results/tables/markers_cluster1.csv, markers_all_clusters.csv,
#         markers_top10_per_cluster.csv
#         results/figures/Violin_plot_for_top_2_markers_of_cluster_1.png,
#         Top_marker_genes.png, Known_marker_genes.png,
#         Top3_markers_DotPlot.png, Top10_markers_heatmap.png
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

cfg <- load_config()
mc <- cfg$markers
log_step("04. Cluster marker genes")

require_files(p_rds(cfg, "03_seurat_clustered.rds"))
obj <- readRDS(p_rds(cfg, "03_seurat_clustered.rds"))
DefaultAssay(obj) <- "RNA"

# ---- Markers of a single cluster (cluster 1) vs all other cells -----------------
if ("1" %in% levels(Idents(obj))) {
  c1 <- FindMarkers(obj, ident.1 = "1", min.pct = mc$min_pct,
                    logfc.threshold = mc$logfc_threshold, verbose = FALSE)
  c1 <- c1[order(c1$p_val_adj, -c1$avg_log2FC), ]
  print(utils::head(c1, 5))
  utils::write.csv(c1, p_table(cfg, "markers_cluster1.csv"))
  save_plot(cfg, VlnPlot(obj, features = rownames(c1)[1:2], pt.size = 0.1),
            "Violin_plot_for_top_2_markers_of_cluster_1", 10, 5)
}

# ---- Markers of every cluster ----------------------------------------------------
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = mc$min_pct,
                          logfc.threshold = mc$logfc_threshold, verbose = FALSE)
markers <- markers[markers$p_val_adj < mc$padj_cutoff, ]
log_msg("Significant positive markers: ", nrow(markers), " (padj < ", mc$padj_cutoff, ")")
utils::write.csv(markers, p_table(cfg, "markers_all_clusters.csv"), row.names = FALSE)

top10 <- as.data.frame(markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
top3 <- as.data.frame(markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top1 <- as.data.frame(markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC))
utils::write.csv(top10, p_table(cfg, "markers_top10_per_cluster.csv"), row.names = FALSE)
print(top3[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], row.names = FALSE)

spectral <- scale_colour_gradientn(colours = rev(brewer.pal(11, "Spectral")))

# Top marker of each cluster on the UMAP.
genes1 <- unique(top1$gene)
p <- FeaturePlot(obj, features = genes1, ncol = 4, order = TRUE)
save_plot(cfg, p, "Top_marker_genes", 16, 4 * ceiling(length(genes1) / 4))

# Canonical markers from the configuration.
known <- unlist(cfg$markers$known)
known <- known[names(known) %in% rownames(obj)]
if (length(known) > 0) {
  plots <- lapply(names(known), function(g) {
    suppressMessages(FeaturePlot(obj, g, order = TRUE) + spectral) +
      ggtitle(paste0(g, ": ", known[[g]]))
  })
  save_plot(cfg, wrap_plots(plots, ncol = 3), "Known_marker_genes", 15,
            4.5 * ceiling(length(plots) / 3))
}

# Dot plot of the top 3 markers per cluster and heatmap of the top 10.
save_plot(cfg, DotPlot(obj, features = unique(top3$gene)) + RotatedAxis() +
            ggtitle("Top 3 markers per cluster"),
          "Top3_markers_DotPlot", max(8, 0.35 * length(unique(top3$gene)) + 3), 6)

hm_assay <- if (identical(cfg$normalization$method, "SCT")) "SCT" else "RNA"
cells <- WhichCells(obj, downsample = 300, seed = cfg$seed)
save_plot(cfg, DoHeatmap(obj, features = unique(top10$gene), cells = cells,
                         assay = hm_assay, size = 3) +
            theme(axis.text.y = element_text(size = 6)),
          "Top10_markers_heatmap", 12, max(8, 0.1 * length(unique(top10$gene)) + 2))

write_session_info(cfg, "04_markers")
