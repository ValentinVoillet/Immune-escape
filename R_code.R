library(tidyverse)
library(here)
library(Seurat)


########################
########################
###---------------------
###--- Patient #1
#- Raw counts
raw_counts <- read.table(file = here("Patient1_raw_countMatrix.csv"), sep = ",", header = TRUE, row.names = 1)

#- Seurat
PBMC <- CreateSeuratObject(counts = raw_counts)
PBMC@meta.data$cellranger.id <- sapply(X = str_split(string = rownames(PBMC@meta.data), pattern = "[.]"), FUN = function(x) x[2])
PBMC@meta.data$time.point <- plyr::mapvalues(x = PBMC@meta.data$cellranger.id, from = c("1", "2"), to = c("2014-02", "2015-06"))
PBMC[["percent.mito"]] <- PercentageFeatureSet(PBMC, pattern = "^MT-")
View(PBMC@meta.data)

#- Normalization and feature selection
PBMC <- NormalizeData(object = PBMC, verbose = FALSE, normalization.method = "LogNormalize")
PBMC <- FindVariableFeatures(object = PBMC, selection.method = "vst", verbose = FALSE)

#- Correction for CDR (nUMI as proxy) and percent.mito
PBMC <- ScaleData(object = PBMC, verbose = TRUE, vars.to.regress = c("nCount_RNA", "percent.mito"), do.par = TRUE, num.cores = 2)

#- PCA
PBMC <- RunPCA(object = PBMC, npcs = 50, verbose = TRUE)
ElbowPlot(object = PBMC, ndims = 50)

#- Clustering
PBMC <- FindNeighbors(object = PBMC, reduction = "pca", dims = 1:20)
PBMC <- FindClusters(object = PBMC, dims.use = 1:20, verbose = TRUE, resolution = .25)

#- Visualization (UMAP)
PBMC <- RunUMAP(object = PBMC, reduction = "pca", dims = 1:20, seed.use = 1216, n.epochs = 400, min.dist = .2, local.connectivity = 20)
DimPlot(object = PBMC, reduction = "umap", group.by = "time.point", pt.size = 0.5)

#- Annotation
PBMC@meta.data$clusters <- PBMC@meta.data$RNA_snn_res.0.25
PBMC@meta.data$annotation.clusters <- plyr::mapvalues(x = PBMC@meta.data$clusters,
                                                      from = c(0:9),
                                                      to = c("AMLs",
                                                             "CD8+ T",
                                                             "CD4+ T",
                                                             "CD14+ Monocytes",
                                                             "CD8+ Transgenic T",
                                                             "NK",
                                                             "B",
                                                             "Proliferating AMLs",
                                                             "PPBP+",
                                                             "CD16+ Monocytes")) %>% as.character
PBMC@meta.data$clusters <- plyr::mapvalues(x = PBMC@meta.data$annotation.clusters,
                                           from = c("AMLs",
                                                    "CD8+ T",
                                                    "CD4+ T",
                                                    "CD14+ Monocytes",
                                                    "CD8+ Transgenic T",
                                                    "NK",
                                                    "B",
                                                    "Proliferating AMLs",
                                                    "PPBP+",
                                                    "CD16+ Monocytes"),
                                           to = paste("Seurat", 1:10))




########################
########################
###---------------------
###--- Patient #2
#- Raw counts
raw_counts <- read.table(file = here("Patient2_raw_countMatrix.csv"), sep = ",", header = TRUE, row.names = 1)

#- Seurat
PBMC <- CreateSeuratObject(counts = raw_counts)
PBMC@meta.data$cellranger.id <- sapply(X = str_split(string = rownames(PBMC@meta.data), pattern = "[.]"), FUN = function(x) x[2])
PBMC@meta.data$time.point <- plyr::mapvalues(x = PBMC@meta.data$cellranger.id, from = c("1", "2"), to = c("2013-12-02", "2013-12-07"))
PBMC[["percent.mito"]] <- PercentageFeatureSet(PBMC, pattern = "^MT-")
View(PBMC@meta.data)

#- Normalization
seurat.list <- SplitObject(PBMC, split.by = "time.point")
for (i in 1:length(seurat.list)) {
  cat("\nNormalization -", i)
  seurat.list[[i]] <- SCTransform(seurat.list[[i]], vars.to.regress = "percent.mito", verbose = TRUE)
}

#- Data Integration
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 2000, verbose = TRUE)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features, verbose = TRUE)
for (i in 1:length(seurat.list)) {
  cat("\nRunning PCA -", i)
  seurat.list[[i]] <- RunPCA(seurat.list[[i]], features = features, verbose = TRUE)
}
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", dims = 1:50, anchor.features = features, verbose = TRUE)
integrated.PBMC <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", dims = 1:50, verbose = TRUE)
DefaultAssay(integrated.PBMC) <- "integrated"

#- Normalization
integrated.PBMC <- NormalizeData(object = integrated.PBMC, assay = "RNA", normalization.method = "LogNormalize", verbose = TRUE)

#- PCA
integrated.PBMC <- RunPCA(object = integrated.PBMC, npcs = 100, verbose = TRUE)
ElbowPlot(integrated.PBMC, ndims = 100)
chosen.nPCs <- 50

#- UMAP
integrated.PBMC <- RunUMAP(object = integrated.PBMC, dims = 1:chosen.nPCs, verbose = TRUE, seed.use = 1234, min.dist = .01)
DimPlot(object = integrated.PBMC, reduction = "umap", group.by = "time.point")

#- Clustering
set.seed(1234)
resolutions <- c(seq(from = .1, to = 1, by = .1), seq(from = 1.5, to = 2, by = .5))
for(i in resolutions) {
  integrated.PBMC <- FindNeighbors(object = integrated.PBMC, dims = 1:chosen.nPCs, verbose = TRUE)
  integrated.PBMC <- FindClusters(object = integrated.PBMC, resolution = i, verbose = TRUE)
}

#- Annotation
integrated.PBMC@meta.data$integrated.clusters <- integrated.PBMC@meta.data$integrated_snn_res.1
integrated.PBMC@meta.data$integrated.annotation.clusters <- plyr::mapvalues(x = integrated.PBMC@meta.data$integrated.clusters,
                                                                            from = c(0:21),
                                                                            to = c("Blasts #1",
                                                                                   "Blasts #2",
                                                                                   "Blasts #3 (MT-genes)",
                                                                                   "CD14+ Monocytes",
                                                                                   "DCs #1",
                                                                                   "Blasts #1",
                                                                                   "CD14+ Monocytes",
                                                                                   "Blasts #1",
                                                                                   "Blasts #1",
                                                                                   "Blasts #1",
                                                                                   "Blasts #2",
                                                                                   "Blasts #1",
                                                                                   "Blasts #1",
                                                                                   "CD14+CD16+ Monocytes",
                                                                                   "T",
                                                                                   "Blasts #1",
                                                                                   "Blasts #1",
                                                                                   "Blasts #3 (MT-genes)",
                                                                                   "Blasts #1",
                                                                                   "NK",
                                                                                   "DCs #2",
                                                                                   "DCs #3")) %>% as.character
integrated.PBMC@meta.data$integrated.clusters <- plyr::mapvalues(x = integrated.PBMC@meta.data$integrated.annotation.clusters,
                                                                      from = c("CD14+CD16+ Monocytes",
                                                                               "CD14+ Monocytes",
                                                                               "DCs #1",
                                                                               "DCs #2",
                                                                               "DCs #3",
                                                                               "Blasts #1",
                                                                               "Blasts #2",
                                                                               "Blasts #3 (MT-genes)",
                                                                               "T",
                                                                               "NK"),
                                                                      to = paste("Seurat", 1:10))

