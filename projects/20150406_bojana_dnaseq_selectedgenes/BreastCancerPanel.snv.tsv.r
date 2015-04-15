library(ggplot2)

setwd("H:/shengquanhu/projects/Jennifer/20150406_bojana_dnaseq_selectedgenes")
inputfiles<-c("BreastCancerPanel.snv.tsv", "BreastCancerPanel_functional.snv.tsv")

for(inputfile in inputfiles){
  data<-read.table(inputfile, sep="\t", header=T)
  
  response<-read.table("response.tsv", sep="\t", header=T, stringsAsFactor=F)
  data$Sample <- factor(data$Sample, levels = response$SampleID)
  
  colors=c("Green", "Red", "Black")
  names(colors)<-c("PCR", "NO-PCR", "NEAR-PCR")
  
  col <- c("red", "black", "blue")
  names(col)<-c("INDEL","SNP","SNP+INDEL")
  
  png(file=paste0(inputfile, ".png"), width=8000, height=16000, res=300)
  print(qplot(Sample, Gene, data=data, fill=SNV, geom="tile") +
          scale_fill_manual(values=c("red", "black", "skyblue")) +
          theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
                axis.text.y = element_text(size=11, face="bold")) + 
          coord_equal())
  dev.off()
  
  mygenes<-read.table("mygenes.txt")
  mydata<-data[data$Gene %in% mygenes$V1,]
  png(file=paste0(inputfile, ".mygene.png"), width=8000, height=2000, res=300)
  print(qplot(Sample, Gene, data=mydata, fill=SNV, geom="tile") +
          scale_fill_manual(values=c("red", "black", "skyblue")) +
          theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold", color=colors[response$Response]),
                axis.text.y = element_text(size=11, face="bold")) + 
          coord_equal())
  dev.off()
}