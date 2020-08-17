library(tidyverse)
library(here)
library(Seurat)

###--- Open Data
raw_counts <- read.csv(file = here("misc", "x090-manuscript", "NCBI", "raw_countMatrix.csv"), row.names = 1) # Download from NCBI
PBMC.x090 <- CreateSeuratObject(counts = raw_counts, min.cells = 0, min.genes = 0) # Quality control has already been done
PBMC.x090[["percent.mito"]] <- PercentageFeatureSet(object = PBMC.x090, pattern = "^MT-")
PBMC.x090@meta.data$cellranger.id <- sapply(X = rownames(PBMC.x090@meta.data), FUN = function(x) str_split(string = x, pattern = "[.]")[[1]][2])
PBMC.x090@meta.data$time.point <- plyr::mapvalues(x = PBMC.x090@meta.data$cellranger.id, from = c("1", "2"), to = c("2014-02", "2015-06"))

###--- Seurat Pipeline
##-- Normalization & Highly Variable Genes
PBMC.x090 <- NormalizeData(object = PBMC.x090, verbose = FALSE, normalization.method = "LogNormalize")
PBMC.x090 <- FindVariableFeatures(object = PBMC.x090, selection.method = "vst", verbose = FALSE)

##-- Correction for CDR (nUMI as proxy) and percent.mito
PBMC.x090 <- ScaleData(object = PBMC.x090, verbose = TRUE, vars.to.regress = c("nCount_RNA", "percent.mito"), do.par = TRUE, num.cores = 2)

##-- PCA
PBMC.x090 <- RunPCA(object = PBMC.x090, npcs = 50, verbose = TRUE)
ElbowPlot(object = PBMC.x090, ndims = 50)
ProjectDim(PBMC.x090, dims.print = 1:20)

##-- Clustering
PBMC.x090 <- FindNeighbors(object = PBMC.x090, reduction = "pca", dims = 1:20)
PBMC.x090 <- FindClusters(object = PBMC.x090, dims.use = 1:20, verbose = TRUE, resolution = .25)

##-- Visualization (UMAP)
PBMC.x090 <- RunUMAP(object = PBMC.x090, reduction = "pca", dims = 1:20, seed.use = 1216, n.epochs = 400, min.dist = .2, local.connectivity = 20)
DimPlot(object = PBMC.x090, reduction = "umap", group.by = "time.point", pt.size = 0.5)
DimPlot(object = PBMC.x090, reduction = "umap", group.by = "RNA_snn_res.0.25", pt.size = 0.5)

##-- Visualization (tSNE)
PBMC.x090 <- RunTSNE(object = PBMC.x090, reduction = "pca", dims = 1:20, seed.use = 1216)
DimPlot(object = PBMC.x090, reduction = "tsne", group.by = "time.point", pt.size = 0.5)
DimPlot(object = PBMC.x090, reduction = "tsne", group.by = "RNA_snn_res.0.25", pt.size = 0.5)

##-- Output (meta.data)
PBMC.x090@meta.data %>% write.csv(file = here("misc", "x090-manuscript", "GitHub", "meta.data.csv"), row.names = TRUE)
