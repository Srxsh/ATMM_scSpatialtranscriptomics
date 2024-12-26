library(dplyr)
library(Seurat)
library(patchwork)
library(openxlsx)
library(ggplot2)

setwd("C:/Users/Srian/OneDrive/Desktop/ATMM_project")
GBM.data <- Read10X(data.dir = "C:/Users/Srian/OneDrive/Desktop/ATMM_project/ATMM_scRNA_data.RData")
load("C:/Users/Srian/OneDrive/Desktop/ATMM_project/ATMM_scRNA_data.RData")

# Initialize the Seurat object with the raw (non-normalized data).
GBM.Vasc.obj <- CreateSeuratObject(counts = GBM.data, project = "GBM_VASC_R", min.cells = 3, min.features = 200)
GBM.Vasc.obj

# QC
GBM_VASC_R[["percent.mt"]] <- PercentageFeatureSet(GBM_VASC_R, pattern = "^MT-")

# Visualize QC metrics as a violin plot
pdf("QC_vlnplot.pdf",width=10,height = 4)
VlnPlot(GBM_VASC_R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

pdf("scatterplot.pdf",width=10,height = 4)
plot1 <- FeatureScatter(GBM_VASC_R, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(GBM_VASC_R, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

GBM.Vasc.obj <- subset(GBM.Vasc.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#NORMALISATION
GBM.Vasc.obj <- NormalizeData(GBM.Vasc.obj, normalization.method = "LogNormalize", scale.factor = 10000)
GBM.Vasc.obj <- NormalizeData(GBM.Vasc.obj)

#highly variable features identification
GBM.Vasc.obj <- FindVariableFeatures(GBM.Vasc.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GBM.Vasc.obj), 10)


# plot variable features with and without labels
pdf("variablefeatures.pdf",width=10,height = 4)
plot1 <- VariableFeaturePlot(GBM.Vasc.obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 / plot2
dev.off()

#SCALING THE DATA
all.genes <- rownames(GBM.Vasc.obj)
GBM.Vasc.obj <- ScaleData(GBM.Vasc.obj, features = all.genes)

#PCA reduction
GBM.Vasc.obj <- RunPCA(GBM.Vasc.obj, features = VariableFeatures(object =GBM.Vasc.obj))
# Examine and visualize PCA results a few different ways
print(GBM.Vasc.obj[["pca"]], dims = 1:5, nfeatures = 5)

pdf("Vizdimloading.pdf",width=10,height=4)
VizDimLoadings(GBM.Vasc.obj, dims = 1:2, reduction = "pca")
dev.off()
pdf("dimplot.pdf",width=10,height=4)
DimPlot(GBM.Vasc.obj, reduction = "pca") + NoLegend()
dev.off()
pdf("dimheatmap.pdf",width=10,height=4)
DimHeatmap(GBM.Vasc.obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(GBM.Vasc.obj, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


#elbowplot
pdf("elbowplot.pdf",width=10,height=4)
ElbowPlot(GBM.Vasc.obj)
dev.off()

#clustering cells
GBM.Vasc.obj <- FindNeighbors(GBM.Vasc.obj, dims = 1:10)
GBM.Vasc.obj <- FindClusters(GBM.Vasc.obj, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(GBM.Vasc.obj), 5)

#UMAP
GBM.Vasc.obj <- RunUMAP(GBM.Vasc.obj, dims = 1:10)

# note that you can set `label = TRUE` or use the LabClusters function to help label
# individual clusters
pdf("umap.pdf",width=10,height=4)
DimPlot(GBM.Vasc.obj, reduction = "umap")
dev.off()

saveRDS(GBM.Vasc.obj, file = "../output/GBM_VASC_R.rds")

#DGE features
# find all markers of cluster 2
cluster2.markers <- FindMarkers(GBM.Vasc.obj, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(GBM.Vasc.obj, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
GBM.Vasc.obj.markers <- FindAllMarkers(GBM.Vasc.obj, only.pos = TRUE)
GBM.Vasc.obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(GBM.Vasc.obj, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

GBM.Vasc.obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
write.xlsx(top10,file = "./allmarkers.xlsx",rowNames=T,colNames=T)

GBM.Vasc.obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write.xlsx(top20,file = "./allmarkers.xlsx",rowNames=T,colNames=T)

pdf("heatmap.pdf",width=10,height=4)
DoHeatmap(GBM.Vasc.obj, features = top10$gene) + NoLegend()
dev.off()

new.cluster.ids <- c("capillary","Putative vSMC", "inflammatoryEC","PC1","Arteriole","TEC type2","artery","aSMC1 &2","Artery&srteriole","vein&venule","PC2","unknwn","SMC unknown","GBM PC")
names(new.cluster.ids) <- levels(GBM.Vasc.obj)
GBM.Vasc.obj <- RenameIdents(GBM.Vasc.obj, new.cluster.ids)

pdf("ANNOTATEDumap.pdf",width=10,height=4)
DimPlot(GBM.Vasc.obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

plot <- DimPlot(GBM.Vasc.obj, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../output/images/GBM.Vasc.obj.jpg", height = 7, width = 12, plot = plot, quality = 50)

saveRDS(GBM.Vasc.obj, file = "../output/GBM_VASC_R.rds")

#remove all other files except GBM.Vasc.obj
rm(list=setdiff(ls(),"GBM.Vasc.obj"))

#correcting batch effects
library(SeuratData)
install.packages("remotes")  # Install remotes if you don't have it
remotes::install_github("satijalab/seurat-data")

# install data set
InstallData("ifnb")



