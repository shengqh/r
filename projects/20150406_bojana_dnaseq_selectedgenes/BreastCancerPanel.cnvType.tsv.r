setwd("H:/shengquanhu/projects/Jennifer/20150406_bojana_dnaseq_selectedgenes")
inputfile<-"BreastCancerPanel.cnv.gene.type.tsv"
data<-read.table(inputfile, sep="\t", header=T)
title<-"Copy Number Variation detected in capture range"

library(reshape)
library(ggplot2)

md<-melt(data,id.vars="Gene")
md$value<-as.character(md$value)
md[is.na(md$value),"value"]<-""
md[md$value == "","value"]<-"NOT_DETECTED"
md$value<-as.factor(md$value)
colnames(md)<-c("Gene", "Sample", "CNV")

response<-read.table("response.tsv", sep="\t", header=T, stringsAsFactor=F)
md$Sample <- factor(md$Sample, levels = response$SampleID)

colors=c("Green", "Red", "Blue")
names(colors)<-c("PCR", "NO-PCR", "NEAR-PCR")

png(file=paste0(inputfile, ".png"), width=6000, height=12000, res=300)
ggplot(md, aes(Sample, Gene))+  
  geom_tile(data=md, aes(fill=CNV), color="white") +
  scale_fill_manual(values=c("red", "blue", "gray")) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold")) +
coord_equal()
dev.off()

mygenes<-read.table("mygenes.txt")
mydata<-md[md$Gene %in% mygenes$V1,]
png(file=paste0(inputfile, ".mygene.png"), width=8000, height=2000, res=300)
ggplot(mydata, aes(Sample, Gene))+  
  geom_tile(data=mydata, aes(fill=CNV), color="white") +
  scale_fill_manual(values=c("red", "gray")) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold")) +
  coord_equal()
dev.off()

