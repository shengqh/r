library("DESeq2")
library("heatmap3")

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

drawHCA<-function(rldselect, labRow=NA, distfun=function(x) as.dist(1 - cor(t(x), use = "pa"))){
    cexCol = max(1.0, 0.2 + 1/log10(ncol(rldselect)))
    cexRow = max(1.0, 0.1 + 1/log10(nrow(rldselect)))
    heatmap3(rldselect, 
             col = hmcols, 
             margins=c(12,5), 
             scale="r", 
             labRow=labRow,
             distfun=distfun,
             cexCol=cexCol,
             cexRow=cexRow)
}

setwd("H:/shengquanhu/projects/Jennifer/20150630_bojana_tnbc/heatmap")
data<-read.table("H:/shengquanhu/projects/Jennifer/20150630_bojana_tnbc/star_genetable/result/20150630_bojana_tnbc_gene.count",sep="\t",header=T, check.names=F, row.names=1)

mydata<-data[,c("3193-BJ-0001","3193-BJ-0004","3193-BJ-0007","3193-BJ-0011")]

designData<-data.frame(Sample=colnames(mydata), Group=1)
dds=DESeqDataSetFromMatrix(countData = mydata,
                           colData = designData,
                           design = ~1)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
assayvsd<-assay(vsd)
colnames(assayvsd)<-colnames(mydata)
rownames(assayvsd)<-rownames(mydata)

gene1<-unlist(read.table("lar1.txt", header=F, stringsAsFactor=F)[,1])
rnames<-rownames(data[data$Feature_gene_name %in% gene1,])
gene1data<-assayvsd[rnames,]
sdv<-apply(gene1data,1,sd)
gene1data<-gene1data[sdv != 0,]
rownames(gene1data)<-data[rownames(gene1data),"Feature_gene_name"]
png(file="lar1.heatmap.png", width=3000, height=3000, res=300)
drawHCA(gene1data)
dev.off()

gene2<-unlist(read.table("lar2.txt", header=F, stringsAsFactor=F)[,1])
rnames<-rownames(data[data$Feature_gene_name %in% gene2,])
gene2data<-assayvsd[rnames,]
sdv<-apply(gene2data,2,sd)
gene2data<-gene2data[sdv != 0,]
rownames(gene2data)<-data[rownames(gene2data),"Feature_gene_name"]
png(file="lar2.heatmap.png", width=3000, height=3000, res=300)
drawHCA(gene2data, labRo=NULL)
dev.off()
