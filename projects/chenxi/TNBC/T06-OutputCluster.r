setwd("/scratch/cqs/shengq1/breastcancer/final/trainingdataset_old_repeat")

load("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.RData")

clst<-res$cluster$cluster

fileconn<-file("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.groups", "w")

for (i in 1:14) {
  #i=1
  pid<-clst[clst==i]
  
  pid<-as.character(names(pid))
  
  writeLines(c(i, pid), fileconn)
}
  
close(fileconn)
  
  

