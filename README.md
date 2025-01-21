# ATMM_scSpatialtranscriptomics
This repository contains the complete codebase for my project, which integrates scRNA-seq data and spatial transcriptomics from human GBM samples. The project investigates spatial heterogeneity in tumor vascular cell types, employing workflows for data preprocessing, integration, differential expression analysis, and cell-cell interaction assays. 
Previously published scRNA seqeuncing data is taken and pre processed (QC, normalisation, variable feature identification). Then, for the purpose of solving for batch effects, integration algorithm is applied. This is then followed by dimensional reduction (using UMAPs for easy visualisation) to obtain the marker gene list and differential gene expression profile.
The spatial data of T cells and Endothelial cell subtypes are obtained from the patient and using likelihood test, they are made comparable. Then using cell cell interaction assays, the spatial architecure and complexity of TME is observed. 


The steps involving scRNA sequencing analysis and visulaisation have been inspired from R toolkit called 'Seurat'

The steps involving spatial transcriptomics analysis is done using codes obtained from another github repository called InSitu Type. 
the order in which the code needs to be followed is as given below:
1. ScRNA_datapreprocessing.R
2. integration_mBE.R
3. DEG_up&downregulated.R
4. DEGcodes2.R
5. GOenrichment.R
6. Insitutype.R
7. spatialdistance.R
   
