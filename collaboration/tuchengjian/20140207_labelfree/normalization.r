library("affy")
library("limma")
library("preprocessCore")
setwd("H:/shengquanhu/projects/tuchengjian/20140207_labelfree")

savedata<-function(data, prefix, type){
  write.csv(data, file=paste0(prefix, type, ".csv"), row.names=T)
  
  png(file=paste0(prefix, type, "_MAplot.png"), width=1000 * ncol(data), height=1000 * ncol(data), res=300)
  mva.pairs(data, ylim=c(-0.5,0.5))
  dev.off()
  
  png(file=paste0(prefix, type, "_boxplot.png"), width=4000, height=4000, res=300)
  boxplot(data)
  dev.off()
}

files<-c("5consecutive_runs.csv", "5discontinues_runs.csv", "Different_load.csv")

file<-files[2]
for(file in files){
  #raw data
  prefix<-gsub("(.*_).+", "\\1", file)
  data<-read.csv(file, header=T, row.names=1)
  data<-as.matrix(log2(data))
  data<-data[,order(colnames(data))]
  savedata(data, prefix, "beforenormalization")
  
  #quantile data
  quan_data<-normalize.quantiles(data, copy=T)
  colnames(quan_data)<-colnames(data)
  rownames(quan_data)<-rownames(data)
  savedata(quan_data, prefix, "quantile")
  
  #total intensity
  sum_total<-apply(data,2,sum)
  sum_total_factor<-mean(sum_total) / sum_total
  sum_total_data<-t(apply(data,1,function(x){
    x * sum_total_factor
  }))
  savedata(sum_total_data, prefix, "total_intensity")
  
  #maximum intensity
  max_total<-apply(data,2,max)
  max_total_factor<-mean(max_total) / max_total
  max_total_data<-t(apply(data,1,function(x){
    x * max_total_factor
  }))
  savedata(max_total_data, prefix, "max_intensity")
  
  #median intensity
  median_total<-apply(data,2,median)
  median_total_factor<-mean(median_total) / median_total
  median_total_data<-t(apply(data,1,function(x){
    x * median_total_factor
  }))
  savedata(median_total_data, prefix, "median_intensity")
  
  #upper-quantile intensity
  upperquantile_total<-apply(data,2,quantile)[4,]
  upperquantile_total_factor<-mean(upperquantile_total) / upperquantile_total
  upperquantile_total_data<-t(apply(data,1,function(x){
    x * upperquantile_total_factor
  }))
  savedata(upperquantile_total_data, prefix, "upperquantile_intensity")
  
  #loess non-linear regression
  loess_data<-normalize.loess(data, log.it=F, maxit=2)
  colnames(loess_data)<-colnames(data)
  rownames(loess_data)<-rownames(data)
  savedata(loess_data, prefix, "loess_ratio")
  
  #limma lowess
  limma_data<-normalizeCyclicLoess(data, iterations=2)
  colnames(limma_data)<-colnames(data)
  rownames(limma_data)<-rownames(data)
  savedata(limma_data, prefix, "limma_lowess")
}
