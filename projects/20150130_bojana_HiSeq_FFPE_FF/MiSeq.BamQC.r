outputdir<-"H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/miseq/QC3bam/result"
inputfile<-"miseq_mapping_quality.txt"
outputfile<-"miseq_mapping_quality.txt.png"

setwd(outputdir)

library(ggplot2)

fp<-read.table(inputfile, header=T, sep="\t", stringsAsFactor=F)
fp<-melt(fp, id="Sample")
colnames(fp)<-c("Sample","Mapped","Reads")

design<-read.table("H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/miseq/miseq_10pairs.txt", sep="\t", header=T, row.names=1)
fp$Patient=design[fp$Sample,"Pair"]
fp$Type=design[fp$Sample,"Sample.Type"]


png(file=outputfile, height=1800, width=1800, res=300)
g<-ggplot(fp, aes(x=Patient, y=Reads, fill=Mapped))+ geom_bar(stat="identity", width=.5)+
  facet_grid(Type ~ .) +
  ggtitle("MiSeq BAM Quality Control Report") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold", color="black"),
        axis.text.y = element_text(size=11, face="bold", color="black"),
        plot.title=element_text(face="bold", size=18))
print(g)
dev.off()

