library(affy)

setwd("/scratch/cqs/shengq1/chenxi/20141101_breastcancer_microarray/trainingset")

celtable<-read.table("Step_07_CelFileList.tsv", sep="\t", header=T, as.is=TRUE)
celtypes = unique(celtable$Type)
celtype<-celtypes[1]
for(celtype in celtypes){
  typedata<-celtable[celtable$Type==celtype,]
  cat("Processing ", celtype, "\n")
  data<-just.rma(filenames=typedata$File, verbose=TRUE, normalize=FALSE)
  rdatafile<-paste0("Step_08_", celtype, "_justRMA.rdata")
  tsvfile<-paste0("Step_08_", celtype, "_justRMA.tsv")
  save(data, file=rdatafile)
  write.exprs(data,file=tsvfile)
  rm(data)
}
