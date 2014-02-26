setwd("H:/shengquanhu/projects/Jennifer/20140218_bojana_MiSeq_HiSeq/genetable/result")
data<-read.table("20140218_bojana_MiSeq_HiSeq_gene.count", header=T, row.names=1, check.names=F)
countdata<-data[,c(2:7)]
x<-1
spcorr<-lapply(c(1:3), function(x){
  cdata<-cbind(countdata[,x], countdata[,x+3])
  cdata<-cdata[rowSums(cdata)>0,]
  cr<-cor.test(cdata[,1], cdata[,2], method="spearman")
  cr$estimate
})

sptable<-data.frame(sample = paste0("Sample", c(1:3)), spcorr = unlist(spcorr))
write.csv(sptable, "Spearman.csv", row.names=F)
