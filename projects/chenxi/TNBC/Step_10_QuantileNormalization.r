library("affy");
library("preprocessCore");

if (.Platform$OS.type == "unix") {
  setwd("/scratch/cqs/shengq1/breastcancer/final/trainingdataset_old_repeat")
}else{
  setwd("H:/shengquanhu/projects/BreastCancer/Dataset")
}

datafile<-"Step_09_expression_commonprobes.tsv"

batch<-read.table("Step_07_CelFileList.tsv", sep="\t", header=T, as.is=T)
cat("reading ", datafile, "...\n")
data<-read.table(datafile,header=T,row.name=1,as.is=T)

#check if the names are matched
myvars <- colnames(data) %in% batch$Name
newdata <- colnames(data)[!myvars]
if(length(newdata) > 0){
  stop("Names between data and batch definition are not matched!")
}

#reorder the samples by dataset to be matched with Step_07_CelFileList.tsv
rdata<-data[,batch$Name]
odata<-2^rdata
#rm(data)
cat("quantile normalization ...\n")
normalize.quantiles.robust(x=odata,copy=F)
qdata<-log2(odata)
rm(odata)

cat("saving ...\n")
normfile<-"Step_10_expression_commonprobes.qnorm"
write.csv(x=qdata,file=normfile)
save(qdata,file=paste0(normfile,".RData"))

rm(qdata)