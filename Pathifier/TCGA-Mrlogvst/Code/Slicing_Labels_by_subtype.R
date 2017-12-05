args <- commandArgs(trailingOnly = TRUE)

###########################################
Path_to_your_Matrix<-c("/home/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/Pathifier/TCGA-Mrlogvst/Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt")
#Path_to_your_Matrix<-args[1]
Path_to_Labels<-c("/home/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/Pathifier/TCGA-Mrlogvst/Data/Labels_Controls_and_Normal_separated_TCGA.txt")
#Path_to_Labels<-args[2]
Path_of_Results<-c("/home/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/Pathifier/TCGA-Mrlogvst/Data/")
#Path_of_Results<-args[3]

Labels<-read.table(Path_to_Labels)

Items<-setdiff(as.character(unique(Labels$Labels)),"Control")
Items<-Items[order(Items)]

List_of_Labels<-list()

for(k in Items){
  Controldf<-data.frame(Labels[which(as.character(Labels$Labels) %in% c("Control")),])
  rownames(Controldf)<-rownames(Labels)[which(as.character(Labels$Labels) %in% c("Control"))]
  colnames(Controldf)<-"Label"
  
  mydf<-data.frame(Labels[which(as.character(Labels$Labels) %in% k),])
  rownames(mydf)<-rownames(Labels)[which(as.character(Labels$Labels) %in% k)]
  colnames(mydf)<-"Label"
  
  List_of_Labels[[k]]<-rbind(Controldf,mydf)
 write.table(List_of_Labels[[k]],file= paste0(Path_of_Results,c("Control_and_"),k,c("_Sliced_Labels.txt")))
}

names(List_of_Labels)<-Items


List_of_Labels[[2]]

BM<-read.table(file="/home/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/Pathifier/TCGA-Mrlogvst/Data/Control_and_ Basal_with_indicator.txt",sep="\t")

BM[1:5,1:5]

sum(colnames(BM)==rownames(List_of_Labels[[1]]))
