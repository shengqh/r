library(affy)

setwd("H:/shengquanhu/projects/BreastCancer/Dataset/")

celtable<-read.table("Step_07_CelFileList.tsv", sep="\t", header=T, as.is=TRUE)
celtypes = unique(celtable$Type)
celtype<-celtypes[1]
for(celtype in celtypes){
  typedata<-celtable[celtable$Type==celtype,]
  datasets<-unique(typedata$Dataset)
  dataset<-datasets[1]
  for(dataset in datasets)
  {
    rdatafile<-paste0("Step_08_Expression_", celtype, "_", dataset, ".rdata")
    cat(rdatafile, "\n")
    if(!file.exists(rdatafile)){
      curdata<-typedata[typedata$Dataset==dataset,]
      data<-ReadAffy(filenames=curdata$File, sampleNames=curdata$Name, verbose=TRUE );
      save(data, file=rdatafile)
    }
  }
}
