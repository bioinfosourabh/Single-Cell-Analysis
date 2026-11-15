#!/usr/bin/env Rscript
# =============================================================================
# 05_annotation_singler.R
#
# Reference-based cell type annotation with SingleR (per cell) using the
# celldex references in annotation$references. The consensus label is the
# pruned label from the first reference, or from the next reference where the
# first one was pruned (low-confidence assignment).
#
# Output: results/rds/05_seurat_annotated.rds
#         results/tables/singler_labels_per_cell.csv,
#         singler_consensus_by_cluster.csv, celltype_counts.csv
#         results/figures/ScoreHeatmap_<reference>.png,
#         <reference>_Annotations.png, Consensus_Annotations.png
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(ggplot2)
})

cfg <- load_config()
log_step("05. Cell type annotation (SingleR)")

require_files(p_rds(cfg, "03_seurat_clustered.rds"))
obj <- readRDS(p_rds(cfg, "03_seurat_clustered.rds"))
DefaultAssay(obj) <- "RNA"
test <- GetAssayData(obj, assay = "RNA", layer = "data")   # log-normalised

label_col <- cfg$annotation$label_column
refs <- reference_names(cfg)
results <- list()
for (r in refs) {
  log_msg("Reference: ", r)
  ref <- load_reference(cfg, r)
  shared <- intersect(rownames(test), rownames(ref))
  log_msg("  shared genes: ", length(shared))
  results[[r]] <- SingleR(test = test, ref = ref, labels = ref[[label_col]])
  obj[[paste0("celltype_", r)]] <- results[[r]]$pruned.labels
  log_msg("  pruned (low-confidence) cells: ", sum(is.na(results[[r]]$pruned.labels)))
}

consensus <- results[[1]]$pruned.labels
for (r in refs[-1]) {
  fill <- is.na(consensus)
  consensus[fill] <- results[[r]]$pruned.labels[fill]
}
obj$celltype_consensus <- consensus

lab_tab <- data.frame(cell = colnames(obj), cluster = Idents(obj),
                      obj@meta.data[, paste0("celltype_", refs), drop = FALSE],
                      celltype_consensus = consensus, check.names = FALSE)
utils::write.csv(lab_tab, p_table(cfg, "singler_labels_per_cell.csv"), row.names = FALSE)
by_cluster <- as.data.frame.matrix(table(cluster = Idents(obj),
                                         label = addNA(factor(consensus))))
utils::write.csv(by_cluster, p_table(cfg, "singler_consensus_by_cluster.csv"))
counts <- as.data.frame(sort(table(consensus, useNA = "ifany"), decreasing = TRUE))
names(counts) <- c("celltype", "cells")
print(counts, row.names = FALSE)
utils::write.csv(counts, p_table(cfg, "celltype_counts.csv"), row.names = FALSE)

# ---- Figures --------------------------------------------------------------------
for (r in refs) {
  save_base_plot(cfg, function() {
    plotScoreHeatmap(results[[r]], main = paste("SingleR scores:", r))
  }, paste0("ScoreHeatmap_", r), 10, 7)
  save_plot(cfg, DimPlot(obj, group.by = paste0("celltype_", r), label = TRUE, repel = TRUE,
                         label.size = 3) + NoLegend() + ggtitle(paste(r, "annotation")),
            paste0(r, "_Annotations"), 8, 6)
}
save_plot(cfg, DimPlot(obj, group.by = "celltype_consensus", label = TRUE, repel = TRUE,
                       label.size = 3) + NoLegend() + ggtitle("Consensus cell type annotation"),
          "Consensus_Annotations", 8, 6)

saveRDS(results, p_rds(cfg, "05_singler_results.rds"))
saveRDS(obj, p_rds(cfg, "05_seurat_annotated.rds"))
write_session_info(cfg, "05_annotation_singler")
