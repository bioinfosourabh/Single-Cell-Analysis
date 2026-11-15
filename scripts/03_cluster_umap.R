#!/usr/bin/env Rscript
# =============================================================================
# 03_cluster_umap.R
#
# Shared-nearest-neighbour graph, Louvain clustering, UMAP and t-SNE.
#
# Output: results/rds/03_seurat_clustered.rds
#         results/tables/cluster_sizes.csv
#         results/figures/UMAP_Clusters.png, t-SNE_Clusters.png
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

cfg <- load_config()
cc <- cfg$clustering
log_step("03. Clustering, UMAP and t-SNE")

require_files(p_rds(cfg, "02_seurat_pca.rds"))
obj <- readRDS(p_rds(cfg, "02_seurat_pca.rds"))
dims <- seq_len(min(cc$dims, ncol(Embeddings(obj, "pca"))))

obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
obj <- FindClusters(obj, resolution = cc$resolution, random.seed = cfg$seed, verbose = FALSE)
obj <- RunUMAP(obj, dims = dims, seed.use = cfg$seed, verbose = FALSE)

sizes <- as.data.frame(table(cluster = Idents(obj)))
print(sizes, row.names = FALSE)
utils::write.csv(sizes, p_table(cfg, "cluster_sizes.csv"), row.names = FALSE)
log_msg(nlevels(Idents(obj)), " clusters (resolution ", cc$resolution, ", PCs 1-", max(dims), ")")

save_plot(cfg, DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
            ggtitle("UMAP of clusters"), "UMAP_Clusters", 8, 6)

if (isTRUE(cc$run_tsne)) {
  perplexity <- min(30, floor((ncol(obj) - 1) / 3))
  obj <- RunTSNE(obj, dims = dims, seed.use = cfg$seed, perplexity = perplexity)
  save_plot(cfg, DimPlot(obj, reduction = "tsne", label = TRUE, repel = TRUE) +
              ggtitle("t-SNE of clusters"), "t-SNE_Clusters", 8, 6)
}

saveRDS(obj, p_rds(cfg, "03_seurat_clustered.rds"))
write_session_info(cfg, "03_cluster_umap")
