library("affy");
library("preprocessCore");
library("genefilter")
library("sva")

setwd("H:/shengquanhu/projects/BreastCancer/Dataset")

#original data
datafile<-"Step_09_expression_commonprobes.tsv";

batch<-read.table("Step_07_CelFileList.tsv", sep="\t", header=T, as.is=T)

#quantile normalization data
qnormfile<-"Step_10_expression_commonprobes.qnorm"
qnormdatafile<-paste0(qnormfile,".RData");
if(!file.exists(qnormdatafile)){
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
  #write.csv(x=qdata,file=normfile)
  save(qdata,file=paste0(normfile,".RData"))
}

#batch normalization data
bnormfile<-"Step_11_expression_commonprobes.bnorm";
bnormdatafile<-paste0(bnormfile,".RData");
if(!file.exists(bnormdatafile)){
  if(!exists("qdata")){
    load(qnormdatafile)
  }

  #batch normalization
  mod<-model.matrix(~1, data=batch)
  
  bdata<-ComBat(dat=qdata, batch=batch$Dataset, mod=mod)
  save(bdata, batch, file=bnormdatafile)

  #write.csv(bdata,file=bnormfile);
  save(bdata,file=bnormdatafile);
  
  rm(qdata)
}

if(!exists("bdata")){
  load(bnormdatafile)
}

#probe names
nm<-rownames(bdata)

#filter gene probe
isGeneProbe<-!(nm %in% grep("^A", nm, value=TRUE))
gpdata<-bdata[isGeneProbe,]
rm(bdata)

#find largest IQR probe for each gene
arrayIQR<-apply(gpdata,1,IQR)
probeName<-rownames(gpdata)
largestProbeName<-findLargest(as.vector(probeName),arrayIQR,"hgu133a")

#present each gene by largest IQR probe
isLargestProbe<-rownames(gpdata) %in% largestProbeName
geneData<-gpdata[isLargestProbe,]
rm(gpdata)
geneNames<-unlist(mget(rownames(geneData),envir=hgu133aSYMBOL))
rownames(geneData)<-geneNames

#save gene data
genefile="Step_12_expression_commonprobes.gene.RData"
save(geneData, file=genefile)

#calculate stdev of each gene
geneDataSd <- apply(geneData, 1, function(x){sd(x)})
write.csv(x=geneDataSd, file="Step_13_expression_commonprobes.gene.sd.csv")

#keep gene with larger stdev
sddata<-geneData[geneDataSd > 0.8,]
rm(geneData)
save(sddata, file="Step_13_expression_commonprobes.gene.sd0.8.RData")

#save gene/array table
write.csv(x=sddata, file="Step_13_expression_commonprobes.gene.sd0.8.csv")

