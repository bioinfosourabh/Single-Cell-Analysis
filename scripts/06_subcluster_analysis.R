#!/usr/bin/env Rscript
# =============================================================================
# 06_subcluster_analysis.R
#
# Two follow-up analyses on annotated sub-populations:
#   A. T/NK cells: markers between the consensus T/NK labels and fine-grained
#      SingleR annotation (label.fine of subclusters$fine_reference).
#   B. Myeloid cells: SCTransform (v2) on the subset and markers between the
#      consensus myeloid labels.
# Sub-populations are selected with the regular expressions in the config.
# Labels with fewer than 3 cells are dropped before marker testing; an
# analysis is skipped if fewer than subclusters$min_cells cells or fewer than
# two labels remain.
#
# Output: results/tables/T_NK_markers.csv, T_NK_fine_annotation_counts.csv,
#         myeloid_markers.csv, myeloid_top10_markers.csv
#         results/figures/T_NK_Cell_Clusters.png, T_NK_Top3_Markers.png,
#         T_NK_FineGrained_Annotations.png, Myeloid_Top_Markers_ViolinPlot.png
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(SingleR)
  library(ggplot2)
})

cfg <- load_config()
sc <- cfg$subclusters
log_step("06. Sub-population analyses (T/NK and myeloid cells)")

require_files(p_rds(cfg, "05_seurat_annotated.rds"))
obj <- readRDS(p_rds(cfg, "05_seurat_annotated.rds"))
DefaultAssay(obj) <- "RNA"
labels <- obj$celltype_consensus

#' Subset cells whose consensus label matches `pattern`, dropping labels
#' with fewer than 3 cells. Returns NULL when the subset is too small.
select_subset <- function(pattern, what) {
  hit <- !is.na(labels) & grepl(pattern, labels, ignore.case = TRUE)
  tab <- table(labels[hit])
  keep_labels <- names(tab)[tab >= 3]
  cells <- colnames(obj)[hit & labels %in% keep_labels]
  log_msg(what, ": ", length(cells), " cells in ", length(keep_labels), " label(s): ",
          paste(sprintf("%s (%d)", names(tab), tab), collapse = ", "))
  if (length(cells) < sc$min_cells || length(keep_labels) < 2) {
    log_msg("  skipped (needs >= ", sc$min_cells, " cells and >= 2 labels)")
    return(NULL)
  }
  s <- subset(obj, cells = cells)
  Idents(s) <- "celltype_consensus"
  s
}

# =============================================================================
# A. T and NK cells
# =============================================================================
t_nk <- select_subset(sc$t_nk_pattern, "T/NK cells")
if (!is.null(t_nk)) {
  mk <- FindAllMarkers(t_nk, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5,
                       verbose = FALSE)
  mk <- mk[mk$p_val_adj < cfg$markers$padj_cutoff, ]
  utils::write.csv(mk, p_table(cfg, "T_NK_markers.csv"), row.names = FALSE)
  print(as.data.frame(mk %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))[, c("cluster", "gene", "avg_log2FC", "p_val_adj")],
        row.names = FALSE)

  save_plot(cfg, DimPlot(t_nk, label = TRUE, repel = TRUE, label.size = 3) + NoLegend() +
              ggtitle("T/NK cells (consensus labels)"), "T_NK_Cell_Clusters", 8, 6)

  top3 <- head(mk[order(-mk$avg_log2FC), "gene"], 3)
  if (length(top3) > 0) {
    save_plot(cfg, FeaturePlot(t_nk, features = top3, ncol = 3, order = TRUE),
              "T_NK_Top3_Markers", 15, 5)
  }

  ref <- load_reference(cfg, sc$fine_reference)
  fine <- SingleR(test = GetAssayData(t_nk, assay = "RNA", layer = "data"),
                  ref = ref, labels = ref$label.fine)
  t_nk$celltype_fine <- fine$pruned.labels
  fine_counts <- as.data.frame(sort(table(fine$pruned.labels, useNA = "ifany"), decreasing = TRUE))
  names(fine_counts) <- c("celltype_fine", "cells")
  print(fine_counts, row.names = FALSE)
  utils::write.csv(fine_counts, p_table(cfg, "T_NK_fine_annotation_counts.csv"), row.names = FALSE)
  save_plot(cfg, DimPlot(t_nk, group.by = "celltype_fine", label = TRUE, repel = TRUE,
                         label.size = 2.5) + NoLegend() +
              ggtitle(paste("Fine-grained T/NK annotation:", sc$fine_reference)),
            "T_NK_FineGrained_Annotations", 9, 7)
  saveRDS(t_nk, p_rds(cfg, "06_seurat_t_nk.rds"))
}

# =============================================================================
# B. Myeloid cells
# =============================================================================
myeloid <- select_subset(sc$myeloid_pattern, "Myeloid cells")
if (!is.null(myeloid)) {
  # SCTransform v2 on the subset; it uses glmGamPoi when installed.
  # An SCT assay inherited from the full object (normalization$method = "SCT")
  # is removed so that the model is fitted on the subset only.
  if ("SCT" %in% SeuratObject::Assays(myeloid)) {
    DefaultAssay(myeloid) <- "RNA"
    myeloid[["SCT"]] <- NULL
  }
  myeloid <- SCTransform(myeloid, vst.flavor = "v2", verbose = FALSE)
  mk <- FindAllMarkers(myeloid, assay = "SCT", only.pos = TRUE, min.pct = 0.5,
                       logfc.threshold = 0.5, verbose = FALSE)
  mk <- mk[mk$p_val_adj < cfg$markers$padj_cutoff, ]
  utils::write.csv(mk, p_table(cfg, "myeloid_markers.csv"), row.names = FALSE)
  top10 <- as.data.frame(mk %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
  utils::write.csv(top10, p_table(cfg, "myeloid_top10_markers.csv"), row.names = FALSE)
  print(as.data.frame(mk %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))[, c("cluster", "gene", "avg_log2FC", "p_val_adj")],
        row.names = FALSE)

  genes <- unique(as.data.frame(mk %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC))$gene)
  genes <- utils::head(unique(c(genes, top10$gene)), 6)
  if (length(genes) > 0) {
    save_plot(cfg, VlnPlot(myeloid, features = genes, assay = "SCT", ncol = 2, pt.size = 0.1),
              "Myeloid_Top_Markers_ViolinPlot", 12, 3.5 * ceiling(length(genes) / 2))
  }
  saveRDS(myeloid, p_rds(cfg, "06_seurat_myeloid.rds"))
}

write_session_info(cfg, "06_subcluster_analysis")
