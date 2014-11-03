library("affy");
library("preprocessCore");

if (.Platform$OS.type == "unix") {
  setwd("/scratch/cqs/shengq1/breastcancer/final/trainingdataset_old_repeat")
}else{
  setwd("H:/shengquanhu/projects/BreastCancer/Dataset")
}

datafile<-"Step_09_expression_commonprobes.tsv"

data<-read.table(datafile,header=T,row.name=1)
tdata<-t(data)
rm(data)

odata<-2^tdata
rm(tdata)
normalize.quantiles.robust(x=odata,copy=F);
qdata<-log2(odata);
rm(odata);

normfile<-paste0(datafile,".qnorm");
save(qdata,file=paste0(normfile,".RData"));

rm(qdata);
