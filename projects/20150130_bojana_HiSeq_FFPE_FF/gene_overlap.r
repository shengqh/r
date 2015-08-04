setwd("H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/")  

tnbcgenes<-read.table("E:/sqh/Dropbox/sciences/projects/20150226_bojana_MiSeq_HiSeq/TNBC_geneid.txt", header=F, stringsAsFactors=F)$V1
hiseq<-unique(read.csv("H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/hiseq/star_deseq2/result/HiSeq_FFPE_VS_FF_DESeq2_sig.csv", header=T, stringsAsFactors=F)$Name)
miseq<-unique(read.csv("H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/miseq/star_deseq2/result/MiSeq_FFPE_VS_FF_DESeq2_sig.csv", header=T, stringsAsFactors=F)$Name)

hi_tnbc<-hiseq[hiseq %in% tnbcgenes]
mi_tnbc<-miseq[miseq %in% tnbcgenes]
hi_mi<-miseq[miseq %in% hiseq]

hi_mi_tnbc<-length(hi_tnbc[hi_tnbc %in% miseq])

hionly<-subset(hiseq, !(hiseq %in% miseq))
hionly<-length(subset(hionly, !(hionly %in% tnbcgenes)))

mionly<-subset(miseq, !(miseq %in% hiseq))
mionly<-length(subset(mionly, !(mionly %in% tnbcgenes)))

tnbconly<-subset(tnbcgenes, !(tnbcgenes %in% hiseq))
tnbconly<-length(subset(tnbconly, !(tnbconly %in% miseq)))

hi_mi_only<-length(hi_mi) - hi_mi_tnbc
hi_tnbc_only<-length(hi_tnbc) - hi_mi_tnbc
mi_tnbc_only<-length(mi_tnbc) - hi_mi_tnbc

library(Vennerable)
png(filename="FFPE_FF_DEGenes_TNBCGenes_overlap.png", width=3000, height =3000, res=300)
common<-Venn(SetNames=c("HiSeq DE genes", "TNBC genes", "MiSeq DE genes"), 
             Weight=c(0, hionly, tnbconly, hi_tnbc_only, mionly, hi_mi_only, mi_tnbc_only, hi_mi_tnbc))
plot(common)
dev.off()
