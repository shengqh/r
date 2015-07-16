setwd("H:/shengquanhu/projects/Jennifer/20150406_bojana_dnaseq_selectedgenes")
inputfile<-"BreastCancerPanel.cnv.gene.tsv"
data<-read.table(inputfile, sep="\t", header=T)
title<-"Number of Copy Number Variation detected in capture range"

library(reshape)

md<-melt(data, id="Gene")
md[is.na(md$value),"value"]<-0
md$Type<-unlist(lapply(md$value, function(x){
  if(x > 0){
    "DUPLICATION"
  }else if(x < 0){
    "DELETION"
  }else{
    ""
  }
}))

md$Type<-factor(md$Type, levels=c("DUPLICATION", "DELETION", ""))
colnames(md)<-c("Gene", "Sample", "Coverage","Type")

response<-read.table("response.tsv", sep="\t", header=T, stringsAsFactor=F)
md$Sample <- factor(md$Sample, levels = response$SampleID)

colors=c("Green", "Red", "Blue")
names(colors)<-c("PCR", "NO-PCR", "NEAR-PCR")

library(ggplot2)

png(file=paste0(inputfile, ".png"), width=6000, height=12000, res=300)
ggplot(md, aes(Sample, Gene))+  
  geom_tile(data=md, aes(fill=Type, alpha=Coverage), color="white") +
  scale_fill_manual(values=c("blue", "red", "white")) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold")) +
  coord_equal()  
dev.off()

mygenes<-read.table("mygenes.txt")
mydata<-md[md$Gene %in% mygenes$V1,]
png(file=paste0(inputfile, ".mygene.png"), width=8000, height=2000, res=300)
ggplot(mydata, aes(Sample, Gene))+  
  geom_tile(data=mydata, aes(fill=Type, alpha=Coverage), color="white") +
  scale_fill_manual(values=c("blue", "white")) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=1, face="bold", color=colors[response$Response]),
        axis.text.y = element_text(size=11, face="bold")) +
  coord_equal()  
dev.off()
