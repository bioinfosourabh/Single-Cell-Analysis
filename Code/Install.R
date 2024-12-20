## Installing Necessary packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install CRAN packages
install.packages('Seurat')
install.packages('RColorBrewer')
install.packages('fastmap')
install.packages('ggplot2')
install.packages('patchwork')
install.packages('magrittr')

# Install Bioconductor packages
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scran")
BiocManager::install("scater")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install(c("BiocGenerics", "SummarizedExperiment"))


## Importing Libraries
library(Seurat)          # For single-cell RNA-seq analysis and visualization  
library(dplyr)           # For data manipulation  
library(patchwork)       # For combining ggplot objects  
library(SingleR)         # For automated cell type annotation  
library(magrittr)        # For the pipe operator  
library(RColorBrewer)    # For color palettes in plots  
library(ggplot2)         # For plot customization  
library(DESeq2)          # For alternative differential expression analysis  
library(scran)           # For robust normalization and marker detection  
library(celldex)         # For reference datasets in cell annotation  
library(SingleCellExperiment) # For managing single-cell data  
library(scater)          # For visualization and quality control  
library(ggrepel)         # For non-overlapping text labels in plots  


