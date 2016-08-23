require(biomaRt)
setwd("H:/shengquanhu/projects/20160705_Evan")
gss<-read.delim("clean_genes.txt", header=F, stringsAsFactors = F)$V1
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
genepos <- lapply(gss, function(x){
  getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
        filters="external_gene_name",
        values = x,
        mart=grch37)})
genepos <- do.call("rbind", genepos)
colnames(genepos)<-c("geneid",	"chr",	"s1",	"s2")
write.table(genepos, file="clean_genes_location.txt", row.names = F, quote = F, sep="\t")

unknowngene<-gss[! (gss %in% genepos$geneid)]
write.table(unknowngene, file="clean_genes_unknown.txt", row.names = F, col.names=F, quote = F, sep="\t")
