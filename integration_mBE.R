library(dplyr)

#Splitting the dataset based on "Patient" metadata
GBM.list <- SplitObject(GBM_VASC_R, split.by = "Patient")

#normalize and identify variable features for each dataset independently
GBM.list <- lapply(X = GBM.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = GBM.list)
GBM.list <- lapply(X = GBM.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Finding the integration anchors
GBM.anchors <- FindIntegrationAnchors(object.list = GBM.list, anchor.features = features, reduction = "rpca")

# Data integration
GBM.combined <- IntegrateData(anchorset = GBM.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(GBM.combined) <- "integrated"

# The standard workflow for visualization and clustering
GBM.combined <- ScaleData(GBM.combined, verbose = FALSE)
GBM.combined <- RunPCA(GBM.combined, npcs = 30, verbose = FALSE)
GBM.combined <- RunUMAP(GBM.combined, reduction = "pca", dims = 1:30)
GBM.combined <- FindNeighbors(GBM.combined, reduction = "pca", dims = 1:30)
GBM.combined <- FindClusters(GBM.combined, resolution = 0.5)
# Visualization
pdf("umap_rerun06.12.pdf", width=10, height=4 )
p1 <- DimPlot(GBM.combined, reduction = "umap")
p1
dev.off()

#Saving the object GBM.combined
saveRDS(GBM.combined, file = "../output/GBM.COMB.rds")

#Default assay is changed to RNA
DefaultAssay(GBM.combined)<-"RNA"

#Marker gene identification with genes being expressed in at least 25% of clusters and minimum fold change of 0.25

GBM.combined.markers <- FindAllMarkers(GBM.combined, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Extracting top 10 and 20 marker genes
GBM.combined.markers%>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

GBM.combined.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
library(openxlsx)
write.xlsx(top20, file = "./top20_markers.xlsx", rowNames = FALSE)

#Heatmap of top 10 marker genes
pdf("heatmap_6.12.24.pdf", width=10,height=4)
DoHeatmap(GBM_VASC_R, features = top10$gene) + NoLegend()
dev.off()

#Annotating clusters
Idents(GBM.integ)<-"seurat_clusters"
new.cluster.ids <- c("Capillary", "Arteriole", "Pericytes","Pericytes", 
                     "Capillary", "Vein& Venule", "Inflammatory EC", "Putative vSMC", 
                     "TEC type2", "aSMC", "Artery", "SMC unknown", "Pericytes","TEC type 1"," GBM PC")

names(new.cluster.ids) <- levels(GBM.integ)
GBM.integ<- RenameIdents(GBM.integ, new.cluster.ids)

#UMAP generation
p1 <- DimPlot(GBM.combined, reduction = "umap")
p1

#Saving the object
saveRDS(GBM.combined, file = "../output/GBM.integ.rds")

#The top marker genes based on the cell type
table(Idents(GBM.combined))
celltype.markers<-FindAllMarkers(GBM.combined,only.pos = T,min.pct = 0.25,logfc.threshold = 0.58)
library(openxlsx)
write.xlsx(celltype.markers, file = "./topmarkers_basedoncelltypes.xlsx", rowNames = TRUE, colnames=TRUE)



# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
pdf("dotplot06.12.24.pdf",width=10,height=4)
DotPlot(GBM.combined, features = c("CAVIN2","SEMA3G","FBLN5","ADGRG6","IL1R1","CX3CL1","ESM1","VWA1","COL3A1","COL18A1","SLC6A13","GPX3","ACTA2","IL6"),group.by = "seurat_clusters") + RotatedAxis()
dev.off()


features <- c("SLC16A1","CAVIN2","EDN1","PLAT","GRASP","TMEM88","IL1R1","VWF","IL6ST","EMP1","TM4SF1","CXCL1","SLC38A5","MYO10","CA2","MGP","S100A6","CLU","ESM1","STC1","COL4A1")


GBM.combined$groups <- sample(c("group1", "group2"), size = ncol(GBM.combined), replace = TRUE)


# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
pdf("Ridgeplot_THY1& COL4A1.pdf",width=15,height=4)
RidgePlot(GBM.combined, features = features, ncol = 2)
dev.off()

# Violin plot - Visualize single cell expression distributions in each cluster
pdf("vlnplt_THY1& COL4A1.pdf",width=10, height=4)
VlnPlot(GBM.combined, features = features)
dev.off()

# Feature plot - visualize feature expression in low-dimensional space
pdf("featureplt_THY1& COL4A1.pdf",width=10, height=4)
FeaturePlot(GBM.combined, features = features)
dev.off()

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
pdf("dotplot_finalcelltype.pdf",width=15, height=4)
DotPlot(GBM.combined, features = features,idents = c("Artery","Arteriole"," Capillary","Venule &Vein","Inflammatory EC","TEC type 1","TEC type 2")) + RotatedAxis()
dev.off()

# Single cell heatmap of feature expression
pdf("heatmap_09.12.pdf", width=10,height=4)
DoHeatmap(subset(GBM.combined, downsample = 100), features = features, size = 3)
dev.off()

#Batch Effect analysis:
pdf("PLOTS(9.12.24)/batcheffect_bfrintegration.pdf", width=8, height=6)
p1 <- DimPlot(GBM_VASC_R, reduction = "umap", group.by = "Patient")
p1
dev.off()