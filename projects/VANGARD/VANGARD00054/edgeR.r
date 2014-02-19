library("preprocessCore");
library("edgeR")

setwd("H:/shengquanhu/projects/VANGARD00054/cuffdiff/result/comparison")

x<-read.table("TREATED_vs_CONTROL.genes.read_group_tracking.count", sep="\t", header=T, row.names=1)
group <- factor(c(rep("CONTROL",3), rep("TREATED",3)))

y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
tb<-et$table
tb$adjustp<-p.adjust(tb$PValue, method="fdr")
tbb<-tb[order(tb$adjustp),]
write.csv(tbb,paste0("TREATED_vs_CONTROL_edgeR_classic.csv"))

sigtbb<-tbb[tbb$adjustp <= 0.2,]

z<-read.table("TREATED_vs_CONTROL.gene_exp.diff", header=T, row.names=1)
fpkm<-read.table("TREATED_vs_CONTROL.genes.read_group_tracking.fpkm",header=T, row.names=1)

names<-rownames(sigtbb)
fpkmsig<-rownames(fpkm) %in% names
fpkmsigtable<-fpkm[fpkmsig,]
fpkmsigtable<-fpkmsigtable[order(rownames(fpkmsigtable)),]

zsig<-rownames(z) %in% names
zsigtable<-z[zsig,]
zsigtable<-zsigtable[order(rownames(zsigtable)),]

rownames(fpkmsigtable)<-zsigtable$gene
write.csv(fpkmsigtable, "TREATED_vs_CONTROL_edgeR_classic.sig.fpkm.csv")

pdf(file="TREATED_vs_CONTROL_edgeR_classic.sig.fpkm.pdf")
hmcols <- colorRampPalette(c("green", "black", "red"))(256)
f<-t(scale(t(fpkmsigtable)))
heatmap(f, col=hmcols,cexCol=0.8)
dev.off()

fpkmsigtable<-read.table("TREATED_vs_CONTROL.gene_exp.diff.sig.fpkm",header=T, row.names=1)
fpkmsigtable<-fpkmsigtable[order(rownames(fpkmsigtable)),]

names<-rownames(fpkmsigtable)

zsig<-rownames(z) %in% names
zsigtable<-z[zsig,]
zsigtable<-zsigtable[order(rownames(zsigtable)),]

rownames(fpkmsigtable)<-zsigtable$gene
write.csv(fpkmsigtable, "TREATED_vs_CONTROL.gene_exp.diff.sig.fpkm.csv")

pdf(file="TREATED_vs_CONTROL.gene_exp.diff.sig.fpkm.pdf")
hmcols <- colorRampPalette(c("green", "black", "red"))(256)
f<-t(scale(t(fpkmsigtable)))
heatmap(f, col=hmcols,cexCol=0.8, cexRow=0.35)
dev.off()
