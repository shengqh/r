library("sfsmisc")
library("genefilter")
library("hgu133a.db")

if (.Platform$OS.type == "unix") {
  setwd("/scratch/cqs/shengq1/breastcancer/final/trainingdataset_old_repeat")
}else{
  setwd("D:/projects/BreastCancer/final");
}

batch<-read.table("breastcancer_affymetrix.tsv.batchdefinition")

dataFile<-"breastcancer_affymetrix.tsv.qnorm";

mdata<-read.csv(dataFile,header=T,row.name=1)

foo<-lm(t(mdata) ~ as.factor(batch$V1))
rm(mdata)

norm_all<-t(foo$res)

rm(foo)

nm <- as.vector(sapply(rownames(norm_all), function(x) substr(x,3,nchar(x)),simplify=TRUE))

rownames(norm_all)<-nm

isGeneProbe<-!(nm %in% grep("^A", nm, value=TRUE))

geneProbeData<-subset(norm_all,isGeneProbe)

rm(norm_all)

arrayIQR<-apply(geneProbeData,1,IQR)

probeName<-rownames(geneProbeData)

largestProbeName<-findLargest(as.vector(probeName),arrayIQR,"hgu133a")

#save largestProbeName
probefile=paste0(dataFile,".bnorm.largestprobe.RData")
save(largestProbeName, file=probefile)

nm<-rownames(geneProbeData)

isLargestProbe<-nm %in% largestProbeName

geneData<-subset(geneProbeData,isLargestProbe)

rm(geneProbeData)

geneNames<-unlist(mget(rownames(geneData),envir=hgu133aSYMBOL))

rownames(geneData)<-geneNames

resultFile=paste0(dataFile,".bnorm.genes")

write.csv(x=geneData,file=resultFile)

resultData=paste0(dataFile,".bnorm.genes.RData")

save(geneData,file=resultData);

rm(geneData)

