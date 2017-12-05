args <- commandArgs(trailingOnly = TRUE)

###########################################
#Path_to_your_Matrix<-c("/home/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/Pathifier/TCGA-Mrlogvst/Data/Control_and_ Basal_with_indicator.txt")
Path_to_your_Matrix<-args[1]
#Path_to_Labels<-c("/home/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/Pathifier/TCGA-Mrlogvst/Data/Control_and_Basal_Sliced_Labels.txt")
Path_to_Labels<-args[2]
#Path_of_Results<-c("/home/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/Pathifier/TCGA-Mrlogvst/Data/")
Path_of_Results<-args[3]
#Tumour_subtype<-"Basal"
Tumour_subtype<-args[4]

#################################################
#################################################
if (!require("vsn")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("vsn", dependencies = TRUE)
  library(vsn)}
if (!require("hexbin")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("hexbin", dependencies = TRUE)
  library(hexbin)}
if (!require("DESeq2")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2", dependencies = TRUE)
  library(DESeq2)}

#################################################

Indi_Matrix<-read.table(Path_to_your_Matrix)
Labels<-read.table(Path_to_Labels)
NORMALS <- Indi_Matrix[1,]
Matrix<-Indi_Matrix[-1,]
Matrix<-round(Matrix)

dds <-DESeqDataSetFromMatrix(countData= Matrix,
colData = Labels,
design= ~Label
)
dds$Label<-relevel(dds$Label,ref = "Control")
colData(dds)


#rld<-rlog(dds,blind=FALSE)
#rldblind<-rlog(dds,blind=TRUE)
vsd <- vst(dds, blind=FALSE)
#vsdblind <- vst(dds, blind=TRUE)


#final_rld<-rbind(NORMALS,assay(rld))
#final_rldblind<-rbind(NORMALS,assay(rldblind))
final_vsd<-rbind(NORMALS,assay(vsd))
#final_vsdblind<-rbind(NORMALS,assay(vsdblind))


#write.table(final_rld,file=paste0(Path_of_Results,c("rld_NOblind"),Tumour_subtype,c(".txt")))
#write.table(final_rldblind,file=paste0(Path_of_Results,c("rld_blind"),Tumour_subtype,c(".txt")))
write.table(final_vsd,file=paste0(Path_of_Results,c("vsd_NOblind"),Tumour_subtype,c(".txt")))
#write.table(final_vsdblind,file=paste0(Path_of_Results,c("vsd_blind"),Tumour_subtype,c(".txt")))
