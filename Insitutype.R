#load required data
cd.subset<-read.xlsx("D:/spatial transcriptomics/NP_avg.xlsx",rowNames=T,colNames=T)
CosMx.Vasc.obj<-readRDS("D:/spatial transcriptomics/CosMx_EC_obj.rds")

#Loading the integrated dataset
load(file = "D:/output/GBM.integ.rds")

#setting identities and subseting data
GBM.integ$subtype<-Idents(GBM.integ)
Idents(GBM.integ)<-"stim"
GBM.integ.f<-subset(GBM.integ,idents=c("Tumour"))
Idents(GBM.integ.f)<-"subtype"

#Using negative probe value as reference
Vasc.neg.avg<-cd.subset[colnames(CosMx_EC_obj),]
rownames(Vasc.neg.avg)<-paste0(Vasc.neg.avg$cell_ID,"_",Vasc.neg.avg$fov,"_",substring(Vasc.neg.avg$sample_id,10))

identical(colnames(CosMx_EC_obj),rownames(Vasc.neg.avg))
[1] TRUE

Vasc.neg<-Vasc.neg.avg$value
names(Vasc.neg)<-rownames(Vasc.neg.avg)
head(Vasc.neg)

#Reference expression profile
Vasc.subcluster.ref.exp<-as.matrix(AverageExpression(GBM.integ.f)[[1]])
dim(Vasc.subcluster.ref.exp)
[1] 27112     2

Vasc.subcluster.ref.exp<-Vasc.subcluster.ref.exp[,-c(3,6,8,10,12)]
################################################################################################################
###################################### Label Transfer with full gene list ######################################
################################################################################################################

vasc.counts <- t(as.data.frame(CosMx_EC_obj@assays$Nanostring@counts))
str(vasc.counts)

#using simulated IF data
immunofluordata <- matrix(rpois(n = nrow(vasc.counts) * 4, lambda = 100), 
                          nrow(vasc.counts))

# perform automatic cohorting:
cohort <- fastCohorting(immunofluordata,
                        gaussian_transform = TRUE) 
# ("Gaussian_transform = TRUE" maps variables to gaussians in order to 
#  place dramatically different variables on the same scale.)
table(cohort)

#Calculating celltype similarity
sup <- insitutypeML(x =vasc.counts,
                    neg=Vasc.neg,
                    cohort = cohort,
                    reference_profiles = Vasc.subcluster.ref.exp)   

#Get results
cluster<-sup[1]
CosMx_EC_obj@meta.data[,"Vasc_sub"]<-cluster[["clust"]]

DimPlot(CosMx_EC_obj,group.by = "Vasc_sub")
Idents(CosMx_EC_obj)<-"Vasc_sub"


Vasc_sub <- CosMx_EC_obj@meta.data[,"Vasc_sub"]
#saving the object
saveRDS(Vasc_sub, file = "Vasc_sub.rds")
