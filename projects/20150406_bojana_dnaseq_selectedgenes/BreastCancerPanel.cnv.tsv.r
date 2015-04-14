setwd("H:/shengquanhu/projects/Jennifer/20150406_bojana_dnaseq_selectedgenes")
data<-read.table("BreastCancerPanel.cnv.gene.tsv", sep="\t", header=T)
outputfile<-"BreastCancerPanel.cnv.gene.png"
title<-"Number of Copy Number Variation detected in capture range"

library(reshape)

md<-melt(data)
md[is.na(md$value),"value"]<-0
#md<-md[!is.na(md$value),]
colnames(md)<-c("Gene", "Sample", "CNVCount")

response<-read.table("response.tsv", sep="\t", header=T, stringsAsFactor=F)
md$Sample <- factor(md$Sample, levels = response$SampleID)

colors=c("Green", "Red", "Blue")
names(colors)<-c("PCR", "NO-PCR", "NEAR-PCR")

library(ggplot2)

png(file=outputfile, width=8000, height=16000, res=300)
ggplot(data=md, aes(x=Sample, y=Gene)) +       
  geom_point(aes(size=CNVCount), colour="darkblue") +
  ggtitle(title) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold"))
dev.off()

mygenes<-read.table("mygenes.txt")
mydata<-md[md$Gene %in% mygenes$V1,]
png(file="mygene.cnv.png", width=8000, height=2000, res=300)
ggplot(data=mydata, aes(x=Sample, y=Gene)) +       
  geom_point(aes(size=CNVCount), colour="darkblue") +
  ggtitle(title) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold"))
dev.off()

