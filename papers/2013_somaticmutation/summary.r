counts<-read.table("H:/shengquanhu/projects/somaticmutation/Figure1.tsv",sep='\t',header=T, row.names=1)
counts<-t(counts)
pdf("H:/shengquanhu/projects/somaticmutation/Figure1.pdf")
par(mar=c(9, 5 ,3 ,2))
barplot(as.matrix(counts), beside=T,  legend = rownames(counts),  col = 1:nrow(counts),
        xlim=c(0, (ncol(counts) + 1) * nrow(counts) + 10),
        ylab = "Somatic Mutation Count",
        args.legend = list(x = (ncol(counts) + 1) * nrow(counts) + 6, y=max(counts) * 1.05, bty = "n"), las = 3)
dev.off()


getintersect<-function(data, names){
  result<-data
  lapply(names, function(name){
    result<<-result[result[,name] != '',]
  })
  result
}

matched<-read.table("H:/shengquanhu/projects/somaticmutation/all_muTect_varscan2_rsmc.site.tsv",sep='\t',header=T)

tcga<-getintersect(matched,"DNA_TCGA")
dna_mutect<-getintersect(matched,"DNA_MUTECT")
dna_varscan2<-getintersect(matched,"DNA_VARSCAN2")
dna_rsmc<-getintersect(matched,"DNA_RSMC")
rna_mutect<-getintersect(matched,"RNA_MUTECT")
rna_varscan2<-getintersect(matched,"RNA_VARSCAN2")
rna_rsmc<-getintersect(matched,"RNA_RSMC")

