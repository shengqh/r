library("affy")
library("limma")
setwd("H:/shengquanhu/projects/tuchengjian/20140123_labelfree_normalization")

savedata<-function(data, prefix, type){
  write.csv(data, file=paste0(prefix, type, ".csv"), row.names=T)
  
  png(file=paste0(prefix, type, "_maplot.png"), width=1000 * ncol(data), height=1000 * ncol(data), res=300)
  mva.pairs(data, ylim=c(-0.5,0.5))
  dev.off()
  
  png(file=paste0(prefix, type, "_boxplot.png"), width=4000, height=4000, res=300)
  boxplot(data)
  dev.off()
}

#raw data
mapdata<-as.matrix(log2(read.csv("MAP_normalization.csv", header=T, row.names=1)))
mapdata<-mapdata[,c(2,3,4,5,1,7,8,9,10,6)]
colnames(mapdata)<-c(paste0("MAP_a", c(1:5)),paste0("MAP_b", c(1:5)))
savedata(mapdata, "MAP_", "beforenormalization")

bsadata<-as.matrix(log2(read.csv("BSA_normalization.csv", header=T, row.names=1)))
colnames(bsadata)<-c(paste0("BSA_", c("025_01", "05_02", "075_03", "1_04",
                                      "025_05", "05_06", "075_07", "1_08",
                                      "025_09", "05_10", "075_11", "1_12",
                                      "025_17", "05_18", "075_19", "1_20")))
bsadata<-bsadata[,order(colnames(bsadata))]
savedata(bsadata, "BSA_", "beforenormalization")
datas<-list()
datas[[1]]<-mapdata
datas[[2]]<-bsadata
prefixes<-c("MAP_", "BSA_")

i<-2
for(i in c(1:2)){
  data<-datas[[i]]
  prefix<-prefixes[i]
  
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
