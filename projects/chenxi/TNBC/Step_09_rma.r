library(affy)

#setwd("H:/shengquanhu/projects/BreastCancer/Dataset/")
setwd("/gpfs21/scratch/cqs/shengq1/chenxi/20141101_breastcancer_microarray/trainingset")

celtable<-read.table("Step_07_CelFileList.tsv", sep="\t", header=T, as.is=TRUE)
celtypes = unique(celtable$Type)
celtype<-celtypes[1]
for(celtype in celtypes){
  typedata<-celtable[celtable$Type==celtype,]
  datasets<-unique(typedata$Dataset)
  dataset<-datasets[1]
  tdata<-NA
  for(dataset in datasets)
  {
    rdatafile<-paste0("Step_8_Expression_", celtype, "_", dataset, ".rdata")
    cat("Reading ", rdatafile, "...\n")
    load(rdatafile)
    if(is.na(tdata)){
      tdata<-data
    }else{
      cat("Merging ", rdatafile, "...\n")
      tdata<-merge(tdata, data)
      rm(data)
    }
  }
  
  cat("Performing rma for ", celtype, "...\n")
  rmadata<-rma(tdata,normalize=FALSE,background=TRUE)
  rm(tdata)
  tsvfile<-paste0("Step_09_Expression_", celtype, ".tsv")
  cat("Saving ", celtype, " to ", tsvfile, "...\n")
  write.exprs(rmadata,file=tsvfile)
  rm(rmadata)
}
