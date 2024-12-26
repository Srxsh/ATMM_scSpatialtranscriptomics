#Differential Gene Expression analysis

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

#Capillary
GBM.combined$celltype.stim <- paste(Idents(GBM.combined), GBM.combined$stim, sep = "_")
GBM.combined$celltype <- Idents(GBM.combined)
Idents(GBM.combined) <- "celltype.stim"
capDEG.identification <- FindMarkers(GBM.combined, ident.1 = "Capillary_Control", ident.2 = "Capillary_Tumour", verbose = FALSE)
head(capDEG.identification, n = 15)
pdf("capillaryDEG.pdf",width=8,height=6)
plots <- VlnPlot(GBM.combined, features = c("INSR","SPP1"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
#Arteriole
arterioleDEG.identification <- FindMarkers(GBM.combined, ident.1 = "Arteriole_Control", ident.2 = "Arteriole_Tumour", verbose = FALSE)
head(arterioleDEG.identification, n = 15)
pdf("arterioleDEG.pdf",width=8,height=6)
plots <- VlnPlot(GBM.combined, features = c("ECSCR","MTRNR2L8"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

#Artery

arteryDEG.identification <- FindMarkers(GBM.combined, ident.1 = "Artery_Control", ident.2 = "Artery_Tumour", verbose = FALSE)
head(arteryDEG.identification, n = 15)
pdf("arteryDEG.pdf",width=8,height=6)
plots <- VlnPlot(GBM.combined, features = c("TIMP1","SPARC"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

#Vein&Venule

VVDEG.identification <- FindMarkers(GBM.combined, ident.1 = "Vein& Venule_Control", ident.2 = "Vein& Venule_Tumour", verbose = FALSE)
head(VVDEG.identification, n = 15)
pdf("VVDEG.pdf",width=8,height=6)
plots <- VlnPlot(GBM.combined, features = c("MGP","ACKR1"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

#Pericytes

pericyteDEG.identification <- FindMarkers(GBM.combined, ident.1 = "Pericytes_Control", ident.2 = " GBM PC_Tumour", verbose = FALSE)
head(pericyteDEG.identification, n = 15)
pdf("pericyteDEG.pdf",width=8,height=6)
plots <- VlnPlot(GBM.combined, features = c("COL3A1","COL4A1"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
#saving files
write.csv(capDEG.identification, "DEG_Capillary.csv", row.names = TRUE)
write.csv(arteryDEG.identification, "DEG_artery.csv", row.names = TRUE)
write.csv(arterioleDEG.identification, "DEG_arteriole.csv", row.names = TRUE)
write.csv(VVDEG.identification, "DEG_VV.csv", row.names = TRUE)
write.csv(pericyteDEG.identification, "DEG_pericyte.csv", row.names = TRUE)

#average gene expression for endothelial and mural cell sybtypes

avg_exp <- AverageExpression(
  object = GBM.integ,
  assays = 'RNA',            # Assay containing the RNA expression data
  features = NULL,           # NULL to include all genes
  return.seurat = FALSE,     # Returns a list instead of a Seurat object
  group.by = "ident",        # Groups by the "ident" column (subtypes)
  layer = "data",            # Use the normalized data layer
  verbose = TRUE             # Prints progress
)
head(avg_exp$RNA)

#saving the files
avg_exp.result<-data.frame(avg_exp)
write.csv(avg_exp.result, "avg_gene_expression.csv", row.names = TRUE)
