library(patchwork)

#Differentially expressed genes in the endothelial and mural cell subtypes
Idents(GBM.combined) <- "cluster.rev"
GBM.combined$celltype.stim <- paste(Idents(GBM.combined), GBM.combined$stim, sep = "_")
GBM.combined$celltype <- Idents(GBM.combined)
Idents(GBM.combined) <- "celltype.stim"
DEG.identification <- FindMarkers(GBM.combined, ident.1 = "Endothelial_Control", ident.2 = "Endothelial_Tumour", verbose = FALSE)
head(muralDEG.identification, n = 15)

#Visualisation
pdf("DEG_endo.pdf", width=8, height=6)
plots <- VlnPlot(GBM.combined, features = c("VWA1","ANGPT2"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
dev.off()

#Comparing the 2 types of endothelial cell subtypes in GBM
tecDEG.identification <- FindMarkers(GBM.integ, ident.1 = "TEC type 1", ident.2 = "TEC type2", verbose = FALSE)
head(tecDEG.identification, n = 15)
write.csv(tecDEG.identification, "DEG_TEC1&2.csv", row.names = TRUE)

# Subset the Seurat object to include only TEC type 1 and TEC type 2
TEC_subset <- subset(GBM.integ, idents = c("TEC type 1", "TEC type2"))

# Define the features of interest
features <- c("ESM1", "STC1","COL4A1","MGP","ANGPT2","MT1E","CAVIN2","SLC39A10","IFI6","MT2A")
features <- c("ESM1", "STC1","PLVAP","COL4A1","MGP","ANGPT2","CYR61","COL4A2","PNP","CD9","MT1E","CAVIN2","SLC38A5","SLC39A10","SLC16A1","MT2A","HOPX","SLC7A5","CSRP2","CD320")

# Generate the violin plot
VlnPlot(TEC_subset, 
        features = features, 
        pt.size = 0) +  
  RotatedAxis()       

#dotplot
pdf("dotplot_TEC1&TEC2_DEG.pdf",width=10,height=6)
DotPlot(GBM.integ, features = features,idents = c("TEC type 1","TEC type2")) + RotatedAxis()
dev.off()
gene.up<-rownames(subset(tecDEG.identification,tecDEG.identification$avg_log2FC>0.58&tecDEG.identification$p_val_adj<0.05))

# Extract downregulated genes (avg_log2FC < -0.58 and adjusted p-value < 0.05)
gene.down <- rownames(subset(tecDEG.identification, tecDEG.identification$avg_log2FC < -0.58 & tecDEG.identification$p_val_adj < 0.05))

# View the downregulated genes
head(gene.down)
head(gene.up)

# Extract upregulated genes (avg_log2FC > 0.58 and adjusted p-value < 0.05)
gene.up <- rownames(subset(tecDEG.identification, tecDEG.identification$avg_log2FC >0.58 & tecDEG.identification$p_val_adj < 0.05))

# View the downregulated and upregulated genes
head(gene.down)
head(gene.up)

