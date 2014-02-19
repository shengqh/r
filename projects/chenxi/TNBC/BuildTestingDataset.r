library("affy");
library("preprocessCore");
library("sfsmisc")
library("genefilter")
library("hgu133a.db")

if (.Platform$OS.type == "unix") {
  setwd("/scratch/cqs/shengq1/breastcancer/final/testingdataset")
  trainingdir<-"/scratch/cqs/shengq1/breastcancer/final/trainingdataset_rdata"
}else{
  setwd("D:/projects/BreastCancer/final");
  trainingdir<-"D:/projects/BreastCancer/final"
}

readcsv<-0;

#original data
datafile<-"breastcancer_test.tsv";

#quantile normalization data
qnormfile<-paste0(datafile,".qnorm");
qnormdatafile<-paste0(qnormfile,".RData");
if(!file.exists(qnormdatafile)){
  #read/transfer data
  data<-read.table(datafile,header=T,row.name=1);
  tdata<-t(data);
  rm(data);
  
  #quantile normalization
  normalize.quantiles.robust(x=tdata,copy=F);
  qnormdata<-log2(tdata);
  rm(tdata);
  
  save(qnormdata,file=qnormdatafile);
  write.csv(qnormdata,file=qnormfile);
}

#batch normalization data
bnormfile<-paste0(qnormfile,".bnorm");
bnormdatafile<-paste0(bnormfile,".RData");
if(!file.exists(bnormdatafile)){
  if(readcsv){
    qnormdata<-read.csv(file=qnormfile,header=T,row.names=1);
  }else{
    load(qnormdatafile);
  }
  
  #batch normalization
  batchdefinitionfile<-paste0(datafile,".batchdefinition");
  batch<-read.table(batchdefinitionfile)
  foo<-lm(t(qnormdata) ~ as.factor(batch$V1))
  rm(qnormdata)
  bnormdata<-t(foo$res)
  rm(foo)
  
  save(bnormdata,file=bnormdatafile);
  write.csv(bnormdata,file=bnormfile);
}

#from probe to genes
genesfile=paste0(bnormfile,".bnorm.genes")
genesdatafile<-paste0(genesfile,".RData");
if(!file.exists(genesdatafile)){
  load(bnormdatafile);
  
  nm <- as.vector(sapply(rownames(bnormdata), function(x) substr(x,3,nchar(x)),simplify=TRUE))
  rownames(bnormdata)<-nm
  
  #filter gene probe
  isGeneProbe<-!(nm %in% grep("^A", nm, value=TRUE))
  geneProbeData<-subset(bnormdata,isGeneProbe)
  rm(bnormdata)

  #load largestProbeName
  load(paste0(trainingdir, "/breastcancer_affymetrix.tsv.qnorm.bnorm.largestprobe.RData"))
  
  #present each gene by largest IQR probe
  nm<-rownames(geneProbeData)
  isLargestProbe<-nm %in% largestProbeName
  geneData<-subset(geneProbeData,isLargestProbe)
  rm(geneProbeData)
  
  geneNames<-unlist(mget(rownames(geneData),envir=hgu133aSYMBOL))
  rownames(geneData)<-geneNames
  
  #save gene data
  write.csv(geneData, file=genesfile);
  save(geneData,file=genesdatafile);
}

#keep only genes whose std >= 0.8
sdfile=paste0(genesfile,".sd0.8");
sddatafile=paste0(sdfile,".RData");
if(!file.exists(sddatafile)){
  load(genesdatafile);
  
  #read gene sd from training dataset
  geneDataSd = read.csv(file=paste0(trainingdir, "/breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd.csv"), row.name=1)
  
  #keep gene with larger stdev
  kdata<-subset(geneData, geneDataSd > 0.8);
  rm(geneData);
  
  write.csv(kdata, file=sdfile);
  save(kdata, file=sddatafile);
}
