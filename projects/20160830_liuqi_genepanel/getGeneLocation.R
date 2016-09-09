setwd("H:/shengquanhu/projects/20160830_liuqi_genepanel/documents")
data<-read.csv("TCPS_Labcorp_genes_20160829.csv", stringsAsFactors = F)

genes<-data$GENE
genes[genes=="RSPO2+FUSION"]="RSPO2"
genes[genes=="RSPO3+FUSION"]="RSPO3"
genes<-c(genes, "EIF3E", "PTPRK")
genes<-genes[order(genes)]

host<-"grch37.ensembl.org"
dataset<-"hsapiens_gene_ensembl"

library(biomaRt)
library(gtools)

mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host=host, dataset=dataset)
allgenepos <- lapply(genes, function(x){
  getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name", "strand"),
        filters="external_gene_name",
        values = x,
        mart=mart)})
allgenepos <- do.call("rbind", allgenepos)
allgenepos$strand<-ifelse(allgenepos$strand==1, "+", "-")
allgenepos$score<-1000
colnames(allgenepos)<-c("chr",	"s1",	"s2", "geneid", "strand", "score")
allgenepos<-allgenepos[,c("chr",	"s1",	"s2", "geneid", "score", "strand")]
write.table(allgenepos, file="TCPS_Labcorp_genes_20160829.all.bed", row.names = F,col.names = F,  quote = F, sep="\t")

genepos<-allgenepos[nchar(allgenepos$chr) < 3,]
genepos<-genepos[order(genepos$chr, genepos$s1 ),]
genepos<-genepos[mixedorder(genepos$chr),]
write.table(genepos, file="TCPS_Labcorp_genes_20160829.bed", row.names = F,col.names = F,  quote = F, sep="\t")

