#Library imports
library(Seurat)
library(dplyr)
library(ggplot2)
library(raster)
library(openxlsx)

#load CosMx vascular and T cell seurat object
CosMx.EC.obj<-readRDS("CosMx_EC_obj.Rds")
CosMx.Tcell.obj<-readRDS("CosMx_Tcell_obj.Rds")

#get endothelial spatial coordinates from SGTU1 patients
STGU1.EC.cor<-list()
for (i in 1:32) {
  fov<-paste0("STGU1_fov",i)
  STGU1.EC.cor[[fov]]<-CosMx_EC_obj@images[[fov]]@boundaries$centroids@coords
  rownames(STGU1.EC.cor[[fov]])<-CosMx_EC_obj@images[[fov]]@boundaries$centroids@cells
}

#get T cell spatial coordinates from SGTU1 patients
STGU1.Tcell.cor<-list()
for (i in 1:32) {
  fov<-paste0("STGU1_fov",i)
  STGU1.Tcell.cor[[fov]]<-CosMx_Tcell_obj@images[[fov]]@boundaries$centroids@coords
  rownames(STGU1.Tcell.cor[[fov]])<-CosMx_Tcell_obj@images[[fov]]@boundaries$centroids@cells
}

#Calculating EC to T cell spatial distance
STGU1_EC_T_distance<-list()

for (i in 1:length(STGU1.EC.cor)) {
  fov<-paste0("STGU1_fov",i)
  for (z in 1:nrow(STGU1.EC.cor[[i]])) {
    
    distance<-as.vector(pointDistance(c(STGU1.EC.cor[[i]][z,1],STGU1.EC.cor[[i]][z,2]), STGU1.Tcell.cor[[i]], lonlat=F))
    STGU1_EC_T_distance[[fov]][[rownames(STGU1.EC.cor[[i]])[z]]]<-distance
    
  }
  
}

#Find endothelials that are spatially close to T cell
STGU1.EC.id.close.to.T<-list()
for (i in 1:length(STGU1.EC.cor)) {
  fov<-paste0("STGU1_fov",i)
  for (z in 1:length(STGU1_EC_T_distance[[i]])) {
    
    STGU1.EC.id.close.to.T[[fov]][[rownames(STGU1.EC.cor[[i]])[z]]]<-min(STGU1_EC_T_distance[[i]][[z]])<150
    
  }
}

STGU1.EC.id.close.to.T.f<-list()
for (i in 1:length(STGU1.EC.id.close.to.T)) {
  
  STGU1.EC.id.close.to.T.f[[i]]<-unlist(STGU1.EC.id.close.to.T[[i]])
  
}

STGU1.EC.id.close.to.T.f<-unlist(STGU1.EC.id.close.to.T.f)
STGU1.EC.id.close.to.T.f<-names(subset(STGU1.EC.id.close.to.T.f,STGU1.EC.id.close.to.T.f=="TRUE"))

##Find endothelials that are spatially far away from T cell
STGU1.EC.id.far.from.T<-setdiff(CosMx_EC_obj@images$Run5753_S2@boundaries$centroids@cells,STGU1.EC.id.close.to.T.f)

#Add EC group information back to seurat object
EC.group<-colnames(CosMx_EC_obj)
EC.group[EC.group%in%c(STGU1.EC.id.close.to.T.f,STGG1.EC.id.close.to.T.f,STGU3.EC.id.close.to.T.f)]<-"EC (Close to T)"
EC.group[EC.group%in%c(STGU1.EC.id.far.from.T,STGG1.EC.id.far.from.T,STGU3.EC.id.far.from.T)]<-"EC (Away from T)"
CosMx_EC_obj@meta.data$EC_group<-EC.group

Idents(CosMx_EC_obj)<-EC.group

#Update EC group information to all cells seurat object
Idents(nano.obj.with.label.transfer, cells = colnames(CosMx_EC_obj)) <- Idents(CosMx_EC_obj)
saveRDS(nano.obj.with.label.transfer, file = "D:/spatial transcriptomics/nano.obj.with.label.transfer.rds")


#VISUALISE 
ImageDimPlot(nano.obj.with.label.transfer, fov = "STGU1_fov1", axes = TRUE, cols = "glasbey")

nano.obj.with.label.transfer$EC_group<-Idents(nano.obj.with.label.transfer)
Idents(nano.obj.with.label.transfer)<-"celltype.new"

Idents(CosMx_EC_obj)<-"subset.obj"

Idents(nano.obj.with.label.transfer, cells = colnames(CosMx_EC_obj)) <- Idents(CosMx_EC_obj)

nano.obj.with.label.transfer$EC_sutype<-Idents(nano.obj.with.label.transfer)
DimPlot(nano.obj.with.label.transfer)

table(nano.obj.with.label.transfer$EC_sutype)

table(nano.obj.with.label.transfer$Sample,nano.obj.with.label.transfer$EC_sutype)

table(nano.obj.with.label.transfer$EC_group,nano.obj.with.label.transfer$EC_sutype)

EC_group<-as.character(as.vector(nano.obj.with.label.transfer$EC_group))
STGU4.EC.id<-nano.obj.with.label.transfer@images$Run5753_S4$centroids@cells
EC_group[EC_group %in% setdiff(EC_group,c("Tumor cell","Astrocyte","Neuron","Oligodendrocyte","Endothelial","Myeloid cell",
                                          "Neutrophil","Mural cell","T cell","B cell","Erythrocyte","Fibroblast",
                                          "EC (Close to T)","EC (Away from T)"))]<-NA


nano.obj.with.label.transfer$EC_group<-EC_group
table(nano.obj.with.label.transfer$EC_group,nano.obj.with.label.transfer$EC_sutype)

#VISUALISE 

ImageDimPlot(nano.obj.with.label.transfer,c("TEC type 1", "TEC type 2", "Tumor cell", "T cell", "B cell"), fov = "STGG1_fov9", axes = TRUE, cols = "glasbey")


# Subset the object to include only the specified cell types
desired_cell_types <- c("TEC type 1", "TEC type2", "Tumor cell", "T cell", "B cell")

# Subset the Seurat object based on the desired cell types
subset_obj <- subset(nano.obj.with.label.transfer, idents = desired_cell_types)

# Plot the selected cell types
ImageDimPlot(subset_obj, fov = "STGG1_fov9", axes = TRUE, cols = "glasbey")

ImageDimPlot(nano.obj.with.label.transfer, fov = "STGG1_fov9", cells = WhichCells(nano.obj.with.label.transfer, idents = c("Tumor cell","T cell", "B cell",
                                                                                   "TEC type 1", "TEC type2")), cols = c("red","blue","lightblue","orange","green"), size = 0.6)


ImageDimPlot(nano.obj.with.label.transfer, fov = "STGG1_fov9", cells = WhichCells(nano.obj.with.label.transfer, idents = c("Tumor cell","T cell", "B cell",
                                                                                                                           "TEC type 1", "TEC type2")), cols = c("red","blue","lightblue","orange","green"), size = 0.6)

library(ggplot2)
library(viridis)
library(RColorBrewer)

# Example: Spatial plot with qualitative colors
pdf("STGG1_fov9_TUMORregion.pdf", width=10, height=4)
ImageDimPlot(
  subset.obj,
  fov = "STGG1_fov14",
  cols = c("Tumor cell" = "red", 
           "TEC type 1" = "blue", 
           "TEC type2" = "green", 
           "B cell" = "orange",
           "T cell"="yellow"),
  axes = TRUE
) 
dev.off()


pdf("STGG1_fov12_TUMORregion.pdf", width=10, height=4)
ImageDimPlot(
  nano.obj.with.label.transfer,
  fov = "STGG1_fov12",
  cols = c("Tumor cell" = "red", 
           "TEC type 1" = "blue", 
           "TEC type2" = "green", 
           "B cell" = "orange",
           "T cell"="yellow"),
  axes = TRUE
) 
dev.off()

pdf("STGG1_fov22_PVniche.pdf", width=10, height=4)
ImageDimPlot(
  nano.obj.with.label.transfer,
  fov = "STGG1_fov22",
  cols = c("Tumor cell" = "red", 
           "TEC type 1" = "blue", 
           "TEC type2" = "green", 
           "B cell" = "orange",
           "T cell"="yellow"),
  axes = TRUE
) 
dev.off()

pdf("STGG1_fov1_TLS.pdf", width=10, height=4)
ImageDimPlot(
  nano.obj.with.label.transfer,
  fov = "STGG1_fov1",
  cols = c("Tumor cell" = "red", 
           "TEC type 1" = "blue", 
           "TEC type2" = "green", 
           "B cell" = "orange",
           "T cell"="yellow"),
  axes = TRUE
) 
dev.off()

pdf("STGG1_fov9_TUMORregion.pdf", width=10, height=4)
ImageDimPlot(
  nano.obj.with.label.transfer,
  fov = "STGG1_fov5",
  cols = c("Tumor cell" = "red", 
           "TEC type 1" = "blue", 
           "TEC type2" = "green", 
           "B cell" = "orange",
           "T cell"="yellow"),
  axes = TRUE
) 
dev.off()

saveRDS(cluster, file = "D:/spatial transcriptomics/cluster.rds")
saveRDS(sup, file = "D:/spatial transcriptomics/sup.rds")



