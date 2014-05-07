library("affy")

savedata<-function(data, prefix, type){
  write.csv(data, file=paste0(prefix, type, ".csv"), row.names=T)
  
  png(file=paste0(prefix, type, "_maplot.png"), width=1000 * ncol(data), height=1000 * ncol(data), res=300)
  mva.pairs(data, ylim=c(-0.5,0.5))
  dev.off()
  
  png(file=paste0(prefix, type, "_boxplot.png"), width=4000, height=4000, res=300)
  boxplot(data)
  dev.off()
}


setwd("H:/shengquanhu/projects/tuchengjian/20140310_labelfree")

#data<-read.table("UPS_BvsC_clearFrame.txt", sep="\t", header=T, row.names=1)
data<-read.table("UPS_6Bvs6C_scaffold_frames.txt", sep="\t", header=T, row.names=1)

#loess non-linear regression
loess_data<-normalize.loess(data, log.it=F, maxit=2)
colnames(loess_data)<-colnames(data)
rownames(loess_data)<-rownames(data)
savedata(loess_data, "UPS_6Bvs6C_scaffold_frames", "_loess_ratio")

