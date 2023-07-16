# Single-Cell-Analysis
# A comprehensive GitHub repository dedicated to single cell data analysis using the R Seurat package


#Installing Necessary packages
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
install.packages('Seurat')
install.packages("RColorBrewer")
install.packages('fastmap')
install.packages('celldex')
browseVignettes("celldex")
BiocManager::install("SingleR")
BiocManager::install("SingleCellExperiment")
BiocManager::install("celldex")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))


library(Seurat)
library(dplyr)
library(patchwork)
library(SingleR)
library(magrittr)
library(RColorBrewer)
library(ggplot2)


#loading dataset
data <- Read10X(data.dir = "C:\\Users\\Sourabh\\OneDrive\\Documents\\R\\single")

#Initializing seurat object 
df <- CreateSeuratObject(counts = data, project = "phd", min.cells = 3, min.features = 200)

data <- NULL

#Quality control metrices
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
df[["percent.rb"]] <- PercentageFeatureSet(df, pattern = "^RP[SL]")

#Visualising QC metrices as violin plot
VlnPlot(df, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1)

#FeatureScatter plots for feature-feature relationships
plot1 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot3 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot3

#filtering based on QC metrices  
df <- subset(df, subset = nCount_RNA < 10000 & nFeature_RNA < 2000 & percent.mt < 10 & percent.rb < 40)
#visualizing filtered data
VlnPlot(df, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1)


#Normalizing the data
df <- NormalizeData(df, normalization.method = "LogNormalize")

#Identification of highly variable features
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)

#Identifying 10 most highly variable genes
top10 <- head(VariableFeatures(df), 10)
top10

#ploting variable features with and without labels
plot1 <- VariableFeaturePlot(df)
plot2 <- LabelPoints(plot = plot1, points = top10, repel =  TRUE)

plot1
plot2

#Scaling the data
all.genes <- rownames(df)
df <- ScaleData(df, features = all.genes)

#linear Dimension reduction
df <- RunPCA(df, features = VariableFeatures(object = df))

#Visualizing PCA results a few different ways
print(df[["pca"]], dims = 1:6, nfeatures = 5)

VizDimLoadings(df, dims = 1:6, nfeatures = 15, reduction = "pca")

DimPlot(df, reduction = "pca")

DimHeatmap(df, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(df)


#Clustering the cells
df <- FindNeighbors(df, dims = 1:10)
df <- FindClusters(df, resolution = 0.5)

#Cluster IDs of the first 5 cells
head(Idents(df), 5)

#Run non-linear dimensional reduction (UMAP)
df <- RunUMAP(df, dims = 1:10, verbose = F)

table(df@meta.data$seurat_clusters)

# individual clusters
DimPlot(df, label = T)

#finding all markers of cluster 1
cluster1.markers <- FindMarkers(df, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5) 

VlnPlot(df, features = c(row.names(cluster1.markers)[1], row.names(cluster1.markers)[2]))

#finding markers for every cluster
df.markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
df.markers

#identifying top 3 marker genes for each cluster
x <- df.markers %>% group_by(cluster) %>% top_n(n=3, wt= avg_log2FC)

#plotting top marker gene of each cluster on UMAP
x <- df.markers %>% group_by(cluster) %>% top_n(n=1, wt= avg_log2FC)
FeaturePlot(df, features = x$gene[1:4])
FeaturePlot(df, features = x$gene[5:8])
FeaturePlot(df, features = x$gene[9:12])

#Plotting marker gene of different cells on UMAP based on previously reported markers and PanglaoDB.
#Plot for Fibroblasts
FeaturePlot(df,"COL6A2") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("COL6A2: Fibroblasts")

#Plot for  T Cells
FeaturePlot(df,"IL7R") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("IL7R: CD4 T Cells")

#Plot for Epithelial Cells
FeaturePlot(df,"KRT14") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("KRT14: Epithelial Cells")


filtered_x <- x[x$cluster == 0, ]
FeaturePlot(df, features = filtered_x$gene[1:3])

filtered_x <- x[x$cluster == 3, ]
FeaturePlot(df, features = filtered_x$gene[1:3])


#Cell type annotation using SingleR
#using Human Primary Cell Atlas Data as reference
ref <- celldex::HumanPrimaryCellAtlasData()
ref
colData(ref)

#converting Seurat object to single cell experiment
sce <- as.SingleCellExperiment(DietSeurat(df))
sce

results.main <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.main)
results.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.fine)

# summary of general cell type annotations
table(results.main$pruned.labels)
table(results.fine$pruned.labels)

#adding annotations to the Seurat object metadata
df@meta.data$results.main <- results.main$pruned.labels
df@meta.data$results.fine <- results.fine$pruned.labels

#visualizing fine-grained annotations
df <- SetIdent(df, value = "results.fine")
DimPlot(df, label = T , repel = T, label.size = 3) + NoLegend()

#visualizing main annotations
df <- SetIdent(df, value = "results.main")
DimPlot(df, label = T , repel = T, label.size = 3) + NoLegend()

levels(df)


#
#identifying the different cell types of T/NK cells
#subseting based on the clustering
subset(x = df, idents = c("T_cells", "NK_cell"), invert = TRUE)

#finding markers for T/NK cells cluster
T_df <- subset(x = df, idents = c("T_cells", "NK_cell"))
T_df.markers <- FindAllMarkers(T_df, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)


#Visualising T/NK cells cluster on UMAP
DimPlot(T_df, label = T)


#3 top genes with highest average expression within T/NK cells.
y <- T_df.markers %>% top_n(n=3, wt= avg_log2FC)
y

#UMAP for top 3 genes with highest average expression within T/NK cells
FeaturePlot(T_df, features = y$gene[1:3])

#Fine-grain Annotation of the cluster to identify different types of cells
sce_T <- as.SingleCellExperiment(DietSeurat(T_df))
results_T.fine <- SingleR(test = sce_T,assay.type.test = 1,ref = ref,labels = ref$label.fine)

#table for Different types of cells in T/NK cells cluster
table(results_T.fine$pruned.labels)
T_NK_cells <- table(results_T.fine$pruned.labels)
write.csv(T_NK_cells, file="T_NK_cells.csv",quote = FALSE,row.names = F)

T_df@meta.data$results_T.fine <- results_T.fine$pruned.labels
T_df <- SetIdent(T_df, value = "results_T.fine")

#Visualizing the T/NK cells 
DimPlot(T_df, label = T , repel = T, label.size = 3) + NoLegend()



##
#Identifying differential marker gene for myeloid cells
##
myeloid <- subset(x = df, idents = c("Monocyte", "Neutrophils", "Macrophage", "Erythroblast", "Platelets"))
myeloid.markers <- FindMarkers(myeloid, ident.1 = "Monocyte", ident.2 = "Neutrophils", ident.3 = "Macrophage", ident.4 = "Erythroblast", ident.5 = "Platelets")
myeloid.markers

DefaultAssay(myeloid) <- "RNA"
myeloid <- NormalizeData(myeloid)
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(myeloid)
myeloid <- ScaleData(myeloid, features = all.genes)

myeloid.markers <- FindAllMarkers(myeloid, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

#differentially marker gene identification for myeloid cells
myeloid_top10_markers.de <- as.data.frame(myeloid.markers %>% top_n(n = 10, wt = avg_log2FC))
myeloid_top10_markers.de

#saving into file
write.csv(myeloid_top10_markers.de, file="myeloid_top10_markers.de.csv",quote = FALSE,row.names = F)

