########################################################################
### This script Get the Differential expresed genes in RNAseq data #####
# Writed by Raúl Alejandro Mejía Pedroza                              ##
########################################################################
#######   Data selected by te user #####################################
########################################################################
exp_mat_path<-c("~/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/DGE/Data/expMatrix_TCGA_cBioPortal_no_males.txt")
Labels_path<-c("~/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/DGE/Data/Labels_Controls_and_Normal_separated_TCGA.txt")
results_path<-c("~/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/DGE/Results-TCGA-Mike/")

mylfctreshold<-1
myp.adj<-0.001
mycoresavaiable<-5
# falta el nombre de las cosas a guardar
# Nota la categoria de referencia se llama "Control" 
################################################
######  Libraries needed    ####################
################################################
if (!require("DESeq2")) {
  biocLite("DESeq2")
  library(DESeq2)}
if (!require("BiocParallel")) {
  biocLite("BiocParallel")
  library(BiocParallel)}


########################################################################
#######   Loading the data #############################################
########################################################################
Exp_Mat <- read.table(exp_mat_path ,header=TRUE , sep= "\t")
Exp_Mat<-as.matrix(Exp_Mat)
AllLabels<-read.table(Labels_path)

#######################################################################
############   Cutting the matrices               #####################
#######################################################################
#  
List_of_matrices<-list()
List_of_dfs <-list()
categories <-unique(as.character(AllLabels$Labels))
categories<-setdiff(categories,"Control")

#let's cut the initial matrix and dataframe
for(k in (categories)){
  # Making the vector to subsetting the data.
  selecteditems<-(grepl("Control",as.character(AllLabels$Labels)) | grepl(k,as.character(AllLabels$Labels)))
  # Subsetting the dataframe "ALLLabels"
  List_of_dfs[[k]]<-data.frame(as.character(AllLabels[selecteditems,]))
  rownames(List_of_dfs[[k]])<-rownames(AllLabels)[selecteditems]
  colnames(List_of_dfs[[k]])<-"condition"
  # # Subsetting the matrix "Exp_Mat"
  List_of_matrices[[k]]<-Exp_Mat[,rownames(AllLabels)[selecteditems]]
}

lapply(List_of_matrices,dim)
table(AllLabels)+112
#######################################################################
############    Building the dds objects          #####################
#######################################################################
dds_list<-list()
for(w in (categories)){
  List_of_dfs[[w]]$condition <- relevel(List_of_dfs[[w]]$condition, ref="Control")
  dds_list[[w]]<- DESeqDataSetFromMatrix(countData = round(List_of_matrices[[w]]),
                                colData = List_of_dfs[[w]],
                                design = ~ condition)  
}

#######################################################################
############   Now the DGE function               #####################
#######################################################################

library("BiocParallel")

time1<-proc.time()
DESeq_list<-list()
results_DESeq_list<-list()
for(x in categories){
  register(MulticoreParam(mycoresavaiable))
  DESeq_list[[x]]<- DESeq(dds_list[[x]],parallel = TRUE)
  register(MulticoreParam(mycoresavaiable))
  results_DESeq_list[[x]]<-results(DESeq_list[[x]],parallel = TRUE)
}
time2 = proc.time()-time1

# Saving the results
save(DESeq_list,file=paste0(results_path,"DESeq_list_TCGA-Mike.RData"))
save(results_DESeq_list,file=paste0(results_path,"results_DESeq_list_TCGA-Mike.RData"))

lfc1_results_DESeq_list<-list()
for(x in categories){
  register(MulticoreParam(mycoresavaiable))
  lfc1_results_DESeq_list[[x]]<-results(DESeq_list[[x]],parallel = TRUE, lfcThreshold= mylfctreshold)
}
time3 = proc.time()-time1

# Saving the results and timming
save(lfc1_results_DESeq_list,file=paste0(results_path,"lfc1_results_DESeq_list_TCGA-Mike.RData"))
write.table(c(time2,time3),file="Timing_runfunctions.txt")


padj10_3_lfc1_results_DESeq_list<-list()
for(k in categories){
  pminus10_3<-grepl("TRUE",lfc1_results_DESeq_list[[k]]$padj < myp.adj)
  padj10_3_lfc1_results_DESeq_list[[k]]<-lfc1_results_DESeq_list[[k]][pminus10_3,]
}
#saving the results
save(padj10_3_lfc1_results_DESeq_list,file=paste0(results_path,"padj10_3_lfc1_results_DESeq_list_TCGA-Mike.RData"))
