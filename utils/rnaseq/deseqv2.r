library("DESeq2")
library("RColorBrewer")
library("gplots")
library("heatmap.plus")

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

setwd("E:/test")
data<-read.table("H:/Shengquanhu/projects/vangard/VANGARD00125_jennifer_rnaseq/htseqcount/result/2059-JP.count", row.names=1, header=T, check.names=F)

hasname <- (! is.numeric(data[1,1]))
if(hasname){
  countData<-data[,c(2:ncol(data))]
}else{
  countData<-data
}

pairs=list("G2_vs_G1" =list("G1"=c("2059-JP-6", "2059-JP-7", "2059-JP-8",  "2059-JP-9"),
                            "G2"=c("2059-JP-3", "2059-JP-4", "2059-JP-5")),
           "G3_vs_G1" =list("G1"=c("2059-JP-6", "2059-JP-7", "2059-JP-8",  "2059-JP-9"),
                            "G3"=c("2059-JP-0", "2059-JP-1", "2059-JP-2")))

pairnames=names(pairs)

for(pairname in pairnames){
  #pairname="G2_vs_G1"
  str(pairname)
  gs=pairs[[pairname]]
  gnames=names(gs)
  g1name=gnames[1]
  g2name=gnames[2]
  g1=gs[[g1name]]
  g2=gs[[g2name]]
  c1=countData[,colnames(countData) %in% g1]
  c2=countData[,colnames(countData) %in% g2]
  pairCountData=cbind(c1, c2)
  
  pairColData=data.frame(condition=factor(c(rep(g1name, ncol(c1)), rep(g2name, ncol(c2)))))
  pairColors<-c(rep("RED", ncol(c1)), rep("BLUE", ncol(c2)))
  
  dds=DESeqDataSetFromMatrix(countData = pairCountData,
                             colData = pairColData,
                             design = ~ condition)
  
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  
  res<-results(dds)
  
  if(hasname){
    res$name<-data[,1]
    res<-res[,c(ncol(res), 1:(ncol(res)-1))]
  }
  
  tbb<-res[order(res$padj),]
  write.csv(as.data.frame(tbb),paste0(pairname, "_DESeq2.csv"))
  
  select<- (!is.na(res$padj)) & (res$padj<0.05) & ((res$log2FoldChange >= 1) | (res$log2FoldChange <= -1))
  
  vsd<-varianceStabilizingTransformation(dds,blind=TRUE)
  vsdmatrix<-as.matrix(assay(vsd))
  vsdselect<-vsdmatrix[select,]
  colnames(vsdselect)<-colnames(pairCountData)
  
  png(filename=paste0(pairname, ".png"), width=4000, height =3000, res=300)
  
  clab<-matrix(c(rep("white", ncol(vsdselect)), pairColors), ncol=2, byrow=FALSE)
  colnames(clab)<-c("", "Group")
  par(mar=c(12, 10, 10, 10))
  heatmap.plus(vsdselect, col = hmcols, ColSideColors = clab, margins=c(10,15))
  
  grid.text(g1name, x = unit(0.05, "npc"), y = unit(0.85, "npc"), just = "left", gp=gpar(fontsize=20, col="RED"))
  grid.text(g2name, x = unit(0.05, "npc"), y = unit(0.80, "npc"), just = "left", gp=gpar(fontsize=20, col="BLUE"))
  
  dev.off()
}




