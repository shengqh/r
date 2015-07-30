setwd("H:/shengquanhu/projects/Jennifer/20150630_bojana_tnbc/association")

library(reshape)
library(ggplot2)

response<-read.table("20150630_bojana_tnbc_45.linear.ped", sep="\t", header=T, stringsAsFactor=F,row.names=1, comment.char = "*")
response<-response[order(response$RESPONSE),]
response$GenoTypeName<-paste0("GenoType_", response$IND_ID)
response$AlleleDepthName<-paste0("AlleleDepth_", response$IND_ID)

inputfile<-"20150630_bojana_tnbc_45.snp.single_bscore.epacts.top5000"
data<-read.table(inputfile, sep="\t", header=T, comment.char = "*")
rownames(data)<-paste0(data[,1], "_", data[,2])

genotypes<-read.table("20150630_bojana_tnbc_45_snp.genotype.tsv", header=T, sep="\t", comment.char = "*", check.names=F)
rownames(genotypes)<-paste0(genotypes[,1],"_",genotypes[,2])
gdata<-genotypes[rownames(data),c(response$GenoTypeName, response$AlleleDepthName)]

fdata<-cbind(data, gdata)

write.csv(fdata, file="bscore_top5000.csv")

fdata$Locus<-rownames(fdata)

gedata<-fdata[c(1:100),c("Locus",response$GenoTypeName)]
colnames(gedata)<-c("Locus",response$IND_ID)

md<-melt(gedata)

colnames(md)<-c("Locus", "Sample", "Genotype")

md$Genotype<-as.factor(md$Genotype)
md$Locus<-factor(md$Locus, levels=fdata$Locus)

colors <- c("PCR"="Green","NO-PCR"="Red","NEAR-PCR"="Blue")

resp<-as.factor(response$RESPONSE)

png(file="20150630_bojana_tnbc_45.snp.single_bscore.epacts.top5000.png", width=4000, height=8000, res=300)
ggplot(md, aes(Sample, Locus))+  
  geom_tile(data=md, aes(fill=Genotype), color="white") +
  scale_fill_manual(values=c("red", "blue", "purple")) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=1, face="bold", color=colors[resp]),
        axis.text.y = element_text(size=11, face="bold")) +
  coord_equal()
dev.off()
