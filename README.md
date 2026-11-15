# Single-cell RNA-seq analysis with Seurat and SingleR

[![pbmc3k-pipeline](https://github.com/bioinfosourabh/single-cell-analysis/actions/workflows/pbmc3k-pipeline.yml/badge.svg)](https://github.com/bioinfosourabh/single-cell-analysis/actions/workflows/pbmc3k-pipeline.yml)

An R workflow for droplet-based (10x Genomics) single-cell RNA-seq data with Seurat (v5) and SingleR. It covers:

1. Quality control: filtering cells on detected genes, UMI counts, and mitochondrial and ribosomal read fractions.
2. Normalisation and selection of highly variable genes (LogNormalize, or SCTransform v2).
3. Dimensional reduction: PCA, UMAP and t-SNE.
4. Graph-based clustering.
5. Cluster marker genes.
6. Marker gene plots.
7. Reference-based cell type annotation with SingleR (Human Primary Cell Atlas and Blueprint/ENCODE).
8. Sub-population analysis of T and NK cells, including fine-grained annotation.
9. Marker genes between myeloid cell types.

## Repository structure

```
single-cell-analysis/
├── config.yaml                  # default: 10x Genomics PBMC 3k
├── config_breast_tumour.yaml    # settings used for the figures in Visualizations/
├── run_pipeline.R               # runs the numbered scripts in order
├── install_packages.R           # installs CRAN and Bioconductor dependencies
├── scripts/
│   ├── utils.R                  # shared helpers (config, logging, figure output)
│   ├── 00_download_data.R       # downloads PBMC 3k
│   ├── 01_qc.R                  # section 1
│   ├── 02_normalize_pca.R       # sections 2-3
│   ├── 03_cluster_umap.R        # section 4
│   ├── 04_markers.R             # sections 5-6
│   ├── 05_annotation_singler.R  # section 7
│   └── 06_subcluster_analysis.R # sections 8-9
├── Visualizations/              # figures shown below (breast tumour data)
└── .github/workflows/           # GitHub Actions: runs the workflow on PBMC 3k
```

Every script reads its parameters (paths, QC thresholds, number of PCs, clustering resolution, references) from a YAML configuration: `config.yaml` by default, or the file given with `--config`. Each script sets the random seed and writes `sessionInfo()` to `results/logs/`.

## Data

**Default dataset: PBMC 3k.** `scripts/00_download_data.R` downloads the 10x Genomics *3k PBMCs from a Healthy Donor* filtered gene-barcode matrices (2,700 cells, Chromium v1, hg19), the dataset used in the Seurat guided clustering tutorial, from [cf.10xgenomics.com](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz). It unpacks them to `data/pbmc3k/filtered_gene_bc_matrices/hg19/`. The QC thresholds in `config.yaml` (200-2,500 genes per cell, < 5% mitochondrial reads) are the ones used for this dataset in the Seurat tutorial.

**Breast tumour data (figures in `Visualizations/`).** The figures shown below were produced from a 10x Genomics count matrix of 14 breast tumour samples (33,694 genes x 44,024 cells). Its barcode prefixes (`sc5rJUQ0xx`) are consistent with the breast cancer data of Qian et al. (2020). These data are not redistributed here: an earlier copy in this repository was a truncated `matrix.mtx` and has been removed. To repeat that analysis, obtain the breast cancer count matrix from the authors' [data portal](https://lambrechtslab.sites.vib.be/en/pan-cancer-blueprint-tumour-microenvironment-0), place `barcodes.tsv`, `genes.tsv` and `matrix.mtx` in `data/breast_tumour/`, and run with `config_breast_tumour.yaml`. That configuration has the QC thresholds used for those figures.

Any 10x directory readable by `Seurat::Read10X()` can be analysed by setting `paths$data_dir` in a copy of `config.yaml`.

## Installation

R (>= 4.3) is required. From the repository root:

```sh
Rscript install_packages.R
```

This installs the missing packages: Seurat, SeuratObject, Matrix, dplyr, ggplot2, patchwork, RColorBrewer and yaml from CRAN; SingleR, celldex, SingleCellExperiment, SummarizedExperiment and glmGamPoi from Bioconductor; and, optionally, presto from GitHub for faster marker tests. The celldex reference datasets are downloaded and cached on first use.

## Usage

```sh
Rscript run_pipeline.R                                     # PBMC 3k: download + steps 01-06
Rscript run_pipeline.R --config config_breast_tumour.yaml  # breast tumour data (see Data)
Rscript run_pipeline.R --from 04 --to 06                   # re-run part of the workflow
```

Outputs go to `results/` (or `paths$results`): `figures/` (PNG), `tables/` (CSV), `rds/` (Seurat objects after each step) and `logs/`. The PBMC 3k run takes about 10-15 minutes on a laptop.

## 1. Quality control (`01_qc.R`)

**Purpose.** Remove low-quality cells before downstream analysis. These include empty or damaged droplets with few detected genes, likely doublets with unusually many genes or UMIs, and cells with a high mitochondrial read fraction, which indicates cytoplasmic RNA loss. Genes detected in fewer than three cells are removed.

```r
data <- Read10X(data.dir = cfg$paths$data_dir)
rownames(data) <- gsub("_", "-", rownames(data))
df <- CreateSeuratObject(counts = data, project = cfg$project, min.cells = 3, min.features = 200)

df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
df[["percent.rb"]] <- PercentageFeatureSet(df, pattern = "^RP[SL]")

VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
        ncol = 4, pt.size = 0.1, layer = "counts")
```

![QC metrics before filtering](Visualizations/Pre_filter.png)

```r
plot1 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot1 + plot2 + plot3
```

![QC metric relationships](Visualizations/QC_Metric_relationships.png)

Cells are then filtered with the thresholds in `qc:` of the configuration. For the breast tumour data these were:

```r
df <- subset(df, subset = nCount_RNA < 75000 & nFeature_RNA < 5000 &
                          percent.mt < 15 & percent.rb < 50)
```

![QC metrics after filtering](Visualizations/After_filtering.png)

## 2. Normalisation and highly variable genes (`02_normalize_pca.R`)

**Purpose.** Remove differences in sequencing depth between cells and select the genes with the highest biological variability, which are used for dimensional reduction.

```r
df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 3000)  # 2000 for PBMC 3k
top10 <- head(VariableFeatures(df), 10)
LabelPoints(plot = VariableFeaturePlot(df), points = top10, repel = TRUE)
```

![Highly variable genes](Visualizations/Highly_Variable_Features.png)

Setting `normalization$method: "SCT"` uses `SCTransform(vst.flavor = "v2")` with the mitochondrial fraction regressed out instead. Log-normalised RNA values are still computed, because marker detection and SingleR use them.

## 3. Principal component analysis (`02_normalize_pca.R`)

**Purpose.** Summarise the scaled expression of the variable genes in a small number of components that capture the main axes of variation between cells.

```r
df <- ScaleData(df, features = rownames(df))
df <- RunPCA(df, features = VariableFeatures(object = df))
VizDimLoadings(df, dims = 1:6, nfeatures = 15, reduction = "pca")
DimPlot(df, reduction = "pca")
DimHeatmap(df, dims = 1:2, cells = 500, balanced = TRUE)
ElbowPlot(df)
```

![PCA loadings](Visualizations/PCA_loadings.png)

![PCA](Visualizations/PCA_Plot.png)

![Heatmap of PC1 and PC2](Visualizations/Dimheatmap_for_first_2_PC.png)

![Elbow plot](Visualizations/ElbowPlot_PCA.png)

The elbow plot guides the number of PCs used for clustering and embedding (`clustering$dims`, 10 in both configurations).

## 4. Clustering, UMAP and t-SNE (`03_cluster_umap.R`)

**Purpose.** Group transcriptionally similar cells and embed them in two dimensions for visualisation. A shared-nearest-neighbour graph is built on the first 10 PCs and partitioned with the Louvain algorithm at resolution 0.5.

```r
df <- FindNeighbors(df, dims = 1:10)
df <- FindClusters(df, resolution = 0.5)
df <- RunUMAP(df, dims = 1:10)
df <- RunTSNE(df, dims = 1:10)
DimPlot(df, reduction = "umap", label = TRUE, repel = TRUE)
```

![UMAP of clusters](Visualizations/UMAP_Clusters.png)

![t-SNE of clusters](Visualizations/t-SNE_Clusters.png)

## 5. Cluster marker genes (`04_markers.R`)

**Purpose.** Identify genes enriched in each cluster relative to all other cells. The test is a Wilcoxon rank-sum test on log-normalised data, and Seurat reports Bonferroni-adjusted p-values. Only positive markers detected in at least 25% of cells of the cluster are kept, with padj < 0.05.

```r
cluster1.markers <- FindMarkers(df, ident.1 = 1, min.pct = 0.25)
VlnPlot(df, features = rownames(cluster1.markers)[1:2])
```

![Top 2 markers of cluster 1](Visualizations/Violin_plot_for_top_2_markers_of_cluster_1.png)

```r
df.markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5) %>%
  filter(p_val_adj < 0.05)
top_markers <- df.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
```

Tables: `markers_cluster1.csv`, `markers_all_clusters.csv`, `markers_top10_per_cluster.csv`.

## 6. Marker gene plots (`04_markers.R`)

**Purpose.** Show where marker genes are expressed in the embedding: the top marker of each cluster, canonical markers of expected cell types (listed under `markers$known` in the configuration), a dot plot of the top three markers per cluster and a heatmap of the top ten.

```r
top_markers_1 <- df.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
FeaturePlot(df, features = top_markers_1$gene)
FeaturePlot(df, "COL6A2") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
```

![Top marker genes](Visualizations/Top_marker_genes.png)

Canonical markers in the breast tumour data: COL6A2 (fibroblasts), IL7R (CD4 T cells) and KRT14 (basal epithelial cells).

![COL6A2](Visualizations/Fibroblasts.png)

![IL7R](Visualizations/T_Cells.png)

![KRT14](Visualizations/Epithelial_Cells.png)

![Top 3 markers of cluster 0](Visualizations/Top_3_marker_genes_for_cluster_0.png)

![Top 3 markers of cluster 3](Visualizations/Top_3_marker_genes_for_cluster_3.png)

For PBMC 3k, the configured markers are IL7R, CCR7, CD8A, MS4A1, CD14, FCGR3A, GNLY, FCER1A and PPBP.

## 7. Cell type annotation with SingleR (`05_annotation_singler.R`)

**Purpose.** Assign each cell the reference cell type whose expression profile it correlates with best, independently of the clustering. Two celldex references are used: the Human Primary Cell Atlas and Blueprint/ENCODE (`label.main`). SingleR prunes low-confidence assignments. The consensus label is the Human Primary Cell Atlas label, or the Blueprint/ENCODE label for cells whose first label was pruned.

```r
ref_human_primary <- celldex::HumanPrimaryCellAtlasData()
ref_blueprint <- celldex::BlueprintEncodeData()
test <- GetAssayData(df, assay = "RNA", layer = "data")   # log-normalised counts

results_primary <- SingleR(test = test, ref = ref_human_primary, labels = ref_human_primary$label.main)
results_blueprint <- SingleR(test = test, ref = ref_blueprint, labels = ref_blueprint$label.main)

results_consensus <- ifelse(!is.na(results_primary$pruned.labels),
                            results_primary$pruned.labels, results_blueprint$pruned.labels)
df$celltype_consensus <- results_consensus
plotScoreHeatmap(results_primary)
```

![SingleR scores, HPCA](Visualizations/ScoreHeatmap_Human_primary_cell_atlas.png)

![SingleR scores, Blueprint/ENCODE](Visualizations/ScoreHeatmap_Blueprint.png)

![Consensus annotation](Visualizations/Consensus_Annotations.png)

![HPCA annotation](Visualizations/HumanPrimaryCellAtlas_Annotations.png)

![Blueprint/ENCODE annotation](Visualizations/BlueprintEncode_Annotations.png)

Neither reference contains tumour-specific epithelial states. In the breast tumour data, labels such as "Gametocytes", "Neurons" or "iPS cells" on epithelial clusters reflect the nearest available reference profile, not the presence of those cell types. Labels should be checked against the cluster markers from sections 5-6. Tables: `singler_labels_per_cell.csv`, `singler_consensus_by_cluster.csv`, `celltype_counts.csv`.

## 8. T and NK cell sub-populations (`06_subcluster_analysis.R`)

**Purpose.** Resolve lymphocyte subsets within the annotated T and NK compartment. Cells with T or NK consensus labels (regular expression `subclusters$t_nk_pattern`) are subset. Markers are computed between their labels (min.pct 0.5, log2FC > 0.5), and the cells are re-annotated with the fine-grained Human Primary Cell Atlas labels (`label.fine`).

```r
T_df <- subset(df, cells = t_nk_cells)
Idents(T_df) <- "celltype_consensus"
T_df.markers <- FindAllMarkers(T_df, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
results_T.fine <- SingleR(test = GetAssayData(T_df, layer = "data"), ref = ref_human_primary,
                          labels = ref_human_primary$label.fine)
```

![T/NK cells](Visualizations/T_NK_Cell_Clusters.png)

![Top 3 T/NK markers](Visualizations/T_NK_Top3_Markers.png)

![Fine-grained T/NK annotation](Visualizations/T_NK_FineGrained_Annotations.png)

Tables: `T_NK_markers.csv`, `T_NK_fine_annotation_counts.csv`.

## 9. Marker genes between myeloid cell types (`06_subcluster_analysis.R`)

**Purpose.** Compare myeloid populations (monocytes, macrophages, dendritic cells, neutrophils; `subclusters$myeloid_pattern`) directly with each other. The subset is re-normalised with SCTransform v2, so that the variance model is fitted on the myeloid cells only, and markers are computed on the SCT assay.

```r
myeloid <- SCTransform(myeloid, vst.flavor = "v2")
myeloid.markers <- FindAllMarkers(myeloid, assay = "SCT", only.pos = TRUE,
                                  min.pct = 0.5, logfc.threshold = 0.5)
myeloid_top10_markers <- myeloid.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
```

![Myeloid markers](Visualizations/Myeloid_Top5_Markers_ViolinPlot.png)

Tables: `myeloid_markers.csv`, `myeloid_top10_markers.csv`. Labels with fewer than three cells are dropped before testing. An analysis is skipped, with a message, if fewer than `subclusters$min_cells` cells or fewer than two labels remain.

## Testing

GitHub Actions (`.github/workflows/pbmc3k-pipeline.yml`) installs the dependencies, downloads PBMC 3k and runs steps 00-06 on every push. The figures and tables are kept as a workflow artifact.

## References

1. Hao Y, Stuart T, Kowalski MH, et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology* 42:293-304. doi:10.1038/s41587-023-01767-y
2. Aran D, Looney AP, Liu L, et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology* 20:163-172. doi:10.1038/s41590-018-0276-y
3. Choudhary S, Satija R (2022). Comparison and evaluation of statistical error models for scRNA-seq. *Genome Biology* 23:27. doi:10.1186/s13059-021-02584-9
4. Qian J, Olbrecht S, Boeckx B, et al. (2020). A pan-cancer blueprint of the heterogeneous tumor microenvironment revealed by single-cell profiling. *Cell Research* 30:745-762. doi:10.1038/s41422-020-0355-0

## Author

Sourabh Kumar · [bioinfosourabh@gmail.com](mailto:bioinfosourabh@gmail.com) · [bioinfosourabh.netlify.app](https://bioinfosourabh.netlify.app) · ORCID [0009-0009-5840-9830](https://orcid.org/0009-0009-5840-9830)

For questions or suggestions, please open an issue or send an email.
