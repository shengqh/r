library("affy");
library("preprocessCore");

if (.Platform$OS.type == "unix") {
  setwd("/scratch/cqs/shengq1/breastcancer/final/trainingdataset_old_repeat")
}else{
  setwd("D:/projects/BreastCancer/final");
}

datafile<-"breastcancer_affymetrix.tsv";

data<-read.table(datafile,header=T,row.name=1);
tdata<-t(data);
rm(data);

normalize.quantiles.robust(x=tdata,copy=F);
ldata<-log2(tdata);
rm(tdata);

normfile<-paste0(datafile,".qnorm");
write.csv(ldata,file=normfile);
save(ldata,file=paste0(normfile,".RData"));

rm(ldata);
