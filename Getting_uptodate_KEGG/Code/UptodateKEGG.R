source("http://bioconductor.org/biocLite.R")
if (!require("gage")) {
  biocLite("gage", ask =FALSE)
  library(gage)
}
kegg_gsets<-kegg.gsets(species = "hsa", id.type = "kegg")
dir.create("../Results")
save(kegg_gsets,file="../Results/uptdateKeggPathways.RData")
