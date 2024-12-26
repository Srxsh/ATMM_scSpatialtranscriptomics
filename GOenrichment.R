#installing required libraries
if(!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

#converting the gene identifiers into Entrez IDs
eg = bitr(gene.up, fromType="SYMBOL",toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg) #up regulated genes
eg = bitr(gene.down, fromType="SYMBOL",toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg) #down regulated genes 
#Gene ontology grouping
ggo <- groupGO(gene     = eg$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

head(ggo)
ggo.result<-data.frame(ggo)

#Gene ontology enrichment analysis 
ego <- enrichGO(gene          = eg$ENTREZID,
                universe      = names(gene.down),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
ego.result<-data.frame(ego)


#saving the results
write.csv(eg, "eg_geneup.csv", row.names = TRUE)
write.csv(ego.result, "egoresult_geneup.csv", row.names = TRUE)

