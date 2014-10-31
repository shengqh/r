setwd("H:/shengquanhu/projects/somaticmutation")

library(ggplot2)
library(reshape2)

data<-read.table("TCGA_muTect_varscan2_rsmc.site.tsv", header=T, sep="\t")

cur<-function(data, scorename, limitnum){
  x<-data[!is.na(data[,scorename]),]
  x<-x[order(x[,scorename], decreasing=T),]
  
  x<-x[c(1:limitnum),]
  count<-nrow(x) / 100
  index<-c(1:100)
  numbers<-round(count*index)
  tcga<-x[x$DNA_TCGA != "",]
  total<-nrow(tcga)
  
  unlist(lapply(index, function(y){
    curnum<-numbers[y]
    xx <- x[c(1:curnum),]
    a<-xx[xx$DNA_TCGA != "",]
    nrow(a) 
    #/ curnum
  }))
}

png("rank_score_DNA.png", width=4000, height=3000, res=300)

rsmc<-cur(data, "DNA_RSMC_score", 5000)
mutect<-cur(data, "DNA_MUTECT_score", 5000)
varscan2<-cur(data, "DNA_VARSCAN2_score", 5000)

df<-data.frame(Index=c(1:100), RSMC=rsmc, MUTECT=mutect, VARSCAN2=varscan2)

data_long <- melt(df, id="Index")

ggplot(data=data_long,
       aes(x=Index, y=value, colour=variable)) +
  geom_line()

dev.off()


png("rank_score_RNA.png", width=4000, height=3000, res=300)
rsmc<-cur(data, "RNA_RSMC_score", 5000)
mutect<-cur(data, "RNA_MUTECT_score", 5000)
varscan2<-cur(data, "RNA_VARSCAN2_score", 5000)

df<-data.frame(Index=c(1:100), RSMC=rsmc, MUTECT=mutect, VARSCAN2=varscan2)

data_long <- melt(df, id="Index")

ggplot(data=data_long,
       aes(x=Index, y=value, colour=variable)) +
  geom_line()
dev.off()