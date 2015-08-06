outputdir<-"H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/miseq/fastqc/result"
inputfile<-"FFPE_FF_MiSeq.FastQC.summary.reads.tsv"
outputfile<-"FFPE_FF_MiSeq.FastQC.summary.reads.tsv.png"

setwd(outputdir)

library(ggplot2)

fp<-read.table(inputfile, header=T, sep="\t", stringsAsFactors=F)

png(file=outputfile, height=1300, width=max(2000, 50 * length(unique(fp$File))), res=300)
g<-ggplot(fp, aes(x=Patient, y=Reads, fill=Type))+ geom_bar(stat="identity", width=.5)+
  facet_grid(Type ~ .) +
  ggtitle("MiSeq FastQC Summary Reads") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold", color="black"),
        axis.text.y = element_text(size=11, face="bold", color="black"),
        plot.title=element_text(face="bold", size=18))
print(g)
dev.off()
