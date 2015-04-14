setwd("H:/shengquanhu/projects/Jennifer/20150406_bojana_dnaseq_selectedgenes")

library(reshape)
library(ggplot2)

genelength<-read.table("BreastCancerPanel.geneLength.tsv", sep="\t", header=T, stringsAsFactor=F, row.names=1)
response<-read.table("response.tsv", sep="\t", header=T, stringsAsFactor=F)
mygenes<-read.table("mygenes.txt")

#######################SNP#########################
inputfile<-"BreastCancerPanel.snp.tsv"
outputfile<-"BreastCancerPanel.snp.png"
mygenefile<-"mygene.snp.png"
title<-"Number of SNP detected per 10,000 bases in capture range"

data<-read.table(inputfile, sep="\t", header=T)
md<-melt(data)

md$value<-md$value / genelength[md$Gene,1] * 10000
colnames(md)<-c("Gene", "Sample", "SNPCount")

md$Sample <- factor(md$Sample, levels = response$SampleID)

colors=c("Green", "Red", "Blue")
names(colors)<-c("PCR", "NO-PCR", "NEAR-PCR")

png(file=outputfile, width=8000, height=16000, res=300)
ggplot(data=md, aes(x=Sample, y=Gene)) +       
  geom_point(aes(size=SNPCount), colour="darkblue") +
  ggtitle(title) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold"))
dev.off()

mydata<-md[md$Gene %in% mygenes$V1,]
png(file=mygenefile, width=8000, height=2000, res=300)
ggplot(data=mydata, aes(x=Sample, y=Gene)) +       
  geom_point(aes(size=SNPCount), colour="darkblue") +
  ggtitle(title) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold"))
dev.off()

#######################InDel#########################
inputfile<-"BreastCancerPanel.indel.tsv"
outputfile<-"BreastCancerPanel.indel.png"
mygenefile<-"mygene.indel.png"
title<-"Number of Insertion/Deletion detected per 10,000 bases in capture range"

data<-read.table(inputfile, sep="\t", header=T)
md<-melt(data)

md$value<-md$value / genelength[md$Gene,1] * 10000
colnames(md)<-c("Gene", "Sample", "SNPCount")

md$Sample <- factor(md$Sample, levels = response$SampleID)

colors=c("Green", "Red", "Blue")
names(colors)<-c("PCR", "NO-PCR", "NEAR-PCR")

png(file=outputfile, width=8000, height=16000, res=300)
ggplot(data=md, aes(x=Sample, y=Gene)) +       
  geom_point(aes(size=SNPCount), colour="darkblue") +
  ggtitle(title) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold"))
dev.off()

mydata<-md[md$Gene %in% mygenes$V1,]
png(file=mygenefile, width=8000, height=2000, res=300)
ggplot(data=mydata, aes(x=Sample, y=Gene)) +       
  geom_point(aes(size=SNPCount), colour="darkblue") +
  ggtitle(title) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold"))
dev.off()
