library(affy)

setwd("H:/shengquanhu/projects/BreastCancer/Dataset/")

celtable<-read.table("Step_07_CelFileList.tsv", sep="\t", header=T, as.is=TRUE)
celtypes = unique(celtable$Type)
celtype<-celtypes[1]
for(celtype in celtypes){
  curdata<-celtable[celtable$Type==celtype,]
  celfiles<-curdata$File
  data<-ReadAffy(filenames=curdata$File, sampleNames=curdata$Name );
  rmadata<-rma(data,normalize=FALSE,background=TRUE);
  rm(data);
  tsvfile<-paste0("Step_8_Expression_", celtype, ".tsv");
  write.exprs(rmadata,file=tsvfile);
  rm(rmadata);
}
