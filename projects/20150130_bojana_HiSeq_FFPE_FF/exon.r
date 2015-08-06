setwd("H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/")  

data<-read.table("H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/hiseq/star_exontable/result/FFPE_FF_HiSeq_exon.count", sep="\t", header=T, row.names=1)

#rownames(data)<-data$Feature
sdata<-data[grepl("ENSG00000175106", data$Feature), c("Feature_gene_name", "IG.064","IG.065")]
mdata<-melt(sdata, id="Feature_gene_name")
ggplot(data=mdata, aes(x=Feature_gene_name, y=value)) + facet_grid(variable~.) +
  geom_bar(stat="identity")

hist(sdata)
