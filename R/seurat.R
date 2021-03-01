###run only one time
install.packages(c("Seurat","cowplot","patchwork","tidyverse","dplyr","RColorBrewer"), dependencies = TRUE)
###

# type ?functionname if you want to know detailed information of each function
# e.g. ?Read10X

### Loading library
library(Seurat)
library(cowplot)
library(patchwork)
library(tidyverse)
library(dplyr)

### this demo data was taken from https://www.nature.com/articles/s41590-019-0312-6

###read 10x raw data(put path of your filtered_feature_bc_matrix)
data <- Read10X(data.dir = "/Volumes/bucket/IshikawaU/Masato/scRNAseq/Gulfiya_RNA/Gulfiya_RNA_Integrated/outs/filtered_feature_bc_matrix")
###initialize data
data <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)

###Quality filtering
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")

VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 10)

#Dimentional reduction
# Run the standard workflow for visualization and clustering
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

Exh.combined <- ScaleData(data, verbose = FALSE)
Exh.combined <- RunPCA(Exh.combined, npcs = 30, verbose = FALSE)
ElbowPlot(Exh.combined,ndims = 30)
# t-SNE and Clustering
Exh.combined <- RunTSNE(Exh.combined, reduction = "pca", dims = 1:20, check_duplicates = FALSE,seed.use = 3)

Exh.combined <- FindNeighbors(Exh.combined, reduction = "pca", dims = 1:16)
# you can change clustering by modifying the resolution
Exh.combined <- FindClusters(Exh.combined,resolution = 0.1)

library(RColorBrewer)

cols <- brewer.pal(5, "Dark2")

p <- DimPlot(Exh.combined, reduction = "tsne", pt.size =0.3,cols = cols)
plot(p)
ggsave("Exhaustion.png", plot=p)

f <- FeaturePlot(Exh.combined, features = c("Pdcd1","Tox", "Tcf7", "Havcr2"), reduction = "tsne",cols = c("lightgrey", "red"), min.cutoff = NA)
ggsave("Exhaustion_feature.png", plot=f)

saveRDS(Exh.combined, file = "Exhaustion.combined.rds")
all.markers <- FindAllMarkers(object = Exh.combined)
write.table(all.markers, "all.markers.csv")

###get_normalized_gene_count_data
Normalized.count <- GetAssayData(Exh.combined, slot = "data")
write.csv(Normalized.count, "LogNormalized.csv")
### normalized and scaled gene count
Zscore <- GetAssayData(Exh.combined, slot ="scale.data")
write.csv(Zscore, "scaled.csv")
### raw count
Raw.count <- GetAssayData(Exh.combined@assays$RNA, slot ="counts")
write.csv(Raw.count, "Raw.csv")

###Vlnplot
VlnPlot(Exh.combined,features = "Pdcd1", slot ="data") 
