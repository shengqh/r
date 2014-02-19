library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")
library("gplots")
library("edgeR")

setwd("H:/shengquanhu/projects/Jennifer/20140206_jennifer_rnaseq_scooter")


myDESeq2<-function(count, group, paired=NULL){
  #      require(DESeq2)
  if (is.null(paired)) {
    dds <- DESeqDataSetFromMatrix(countData = count,
                                  colData = data.frame(group=as.factor(group)),
                                  design = formula(~group))
  } else {
    dds <- DESeqDataSetFromMatrix(countData = count,
                                  colData = data.frame(paired=as.factor(paired), group=as.factor(group)),
                                  design = ~paired+ group)
  }
  dds <- DESeq(dds)
  dds <- results(dds)
  return(dds)
}

myedgeR2<-function(count, group,paired=NULL){
  #      require(edgeR)
  y<-DGEList(count)
  group<-as.factor(group)
  if (is.null(paired)) {
    design <- model.matrix(~group)
  } else {
    paired<-as.factor(paired)
    design <- model.matrix(~paired+group)
  }
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTrendedDisp(y,design)
  y <- estimateGLMTagwiseDisp(y,design)
  fit<-glmFit(y, design)
  lrt<-glmLRT(fit, coef=c(-1,1,0,0,0))
  
  edgeR <- topTags(lrt ,nrow(count) )
  head(edgeR$table)
  
  count["ENSG00000115616",]
  
  return(edgeR$table)
}


files<-data.frame(software=c("tophat", "star"), file=c("2059-JP_gene_count_table_all_tophat.txt", "2059-JP_gene_count_table_all_star.txt"), stringsAsFactors=F)

colData<-read.csv("BRE-0904_test_metadata.csv")

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

i<-1
for(i in c(1:nrow(files))){
  software<-files$software[i]
  file<-files$file[i]
  countData<-read.table(file, row.names=1, header=T, check.names=F)
  colnames(countData)<-apply(colData, 1, function(x){
    paste0(substr(x["sequencing_id"], 6, 10), "_P", substr(x["patient_id"], 12, 13))
  }) 
  
  #different expression analysis
  
  dds=DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~1)
  
  colnames(dds)<-colnames(countData)
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  write.csv(assay(vsd), file=paste0("2059-JP-",software, "-DESeq2-vsd.csv"))
  
  assayvsd<-assay(vsd)
  vsdiqr<-apply(assayvsd, 1, IQR)
  assayvsd<-assayvsd[order(vsdiqr, decreasing=T),]
  
  pairColors<-rep("red", 10)
  pairColors[colData$tissue_processing == "whole"]<-"blue"
  
  png(file=paste0("2059-JP-",software,"-heatmap-top200iqr.png"), width=4000, height=4000, res=300)
  heatmap3(assayvsd[1:200,], col = hmcols, ColSideColors = pairColors, ColSideLabs="Processing", margins=c(12,5), scale="r", dist=dist, labRow="",
           legendfun=function() showLegend(legend=c("LCM", "whole"),col=c("red","blue"),cex=1.5,x="center"))
  dev.off()
  
  png(file=paste0("2059-JP-",software,"-sampledist.png"), width=4000, height=4000, res=300)
  distsRL <- dist(t(assay(vsd)))
  heatmap.2(as.matrix(distsRL), trace="none", col = rev(hmcols), margin=c(13, 13))
  dev.off()
  
  indecies<-c(2,3,5:10)
  dataProcess<-countData[,indecies]
  dataProcess<-dataProcess[apply(dataProcess, 1, max) > 0,]
    
  dataProcessColors<-colData[indecies,]
  dataProcessColors$patient_id<-as.factor(as.character(dataProcessColors$patient_id))
  ddsProcess<-DESeqDataSetFromMatrix(countData = dataProcess,
                                      colData = dataProcessColors,
                                      design = ~ patient_id + tissue_processing)
  ddsProcess<-DESeq(ddsProcess)
  
#   png(file=paste0("2059-JP-",software,"-WholeVsLCM-MAplot.png"), width=4000, height=4000, res=300)
#   plotMA(ddsProcess,ylim=c(-2,2),main="DESeq2")
#   dev.off()
  
  res<-results(ddsProcess)
  res<-res[order(res$padj),]
  write.csv(res, paste0("2059-JP-",software,"-WholeVsLCM-DE.csv"))
  
  indecies<-c(1:6)
  dataProcess<-countData[,indecies]
  dataProcess<-dataProcess[apply(dataProcess, 1, max) > 0,]

  dataProcessColors<-colData[indecies,]
  dataProcessColors$patient_id<-as.factor(as.character(dataProcessColors$patient_id))
  ddsProcess<-DESeqDataSetFromMatrix(countData = dataProcess,
                                     colData = dataProcessColors,
                                     design = ~ tnbc_subtype)
  ddsProcess<-DESeq(ddsProcess)
  
#   png(file=paste0("2059-JP-",software,"-SubType-MAplot.png"), width=4000, height=4000, res=300)
#   plotMA(ddsProcess,ylim=c(-2,2),main="DESeq2")
#   dev.off()
  
  res<-results(ddsProcess)
  res<-res[order(res$padj),]
  write.csv(res, paste0("2059-JP-",software,"-SubType-DE.csv"))
}

# count<-dataProcess
# group<-c(1,1,1,1,0,0,0,0)
# paired<-c(1,2,3,4,1,2,3,4)
# 
# Subject<-dataProcessColors$patient_id
# Treat<-dataProcessColors$tissue_processing
# design <- model.matrix(~Subject+Treat)
# 
# ddsDeseq2<-myDESeq2(dataProcess, group, paired)
# ddsEdgeR<-myedgeR2(dataProcess, group, paired)
