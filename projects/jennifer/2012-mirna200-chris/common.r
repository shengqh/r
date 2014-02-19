library(stats)
library(graphics)
library(reshape)
library(edgeR)
library(DESeq)
library(heatmap.plus)
library(gplots)
library(xtable)
library(Heatplus)

setwd("H:/shengquanhu/projects/miRNA")

#load csv file data
loadData<-function(filename){
  read.csv(filename, header=T, row.names=1, check.names=F)
}

#save csv file data
saveData<-function(x, filename){
  write.csv(x, file=filename)
}

#get tumor that each sample belongs to
getSampleTumors<-function(x){
  out<-matrix ( unlist( strsplit(colnames(x), "_") ), byrow=T, ncol=2)
  return(out[ ,1])
}

#get all tumor names
getTumorNames<-function(x){
  samplenames<-getSampleTumors(x)
  nametable<-table(samplenames)
  result<-names(nametable)
  return (result)
}

#build color table of cancers
getTumorColors<-function(x){
  tumornames<-getTumorNames(x)
  colors<-rainbow(length(tumornames))
  names(colors)<-tumornames  
  return (colors)
}

#build color of each sample based on groupname
getSampleColors<-function(x){
  tumorColors<-getTumorColors(x)
  tumorGroups<-getSampleTumors(x)
  return (tumorColors[tumorGroups])
}

#get file name based on sample type and value type
getFile<-function(stype, vtype, extension = ""){
  name<-paste0("TCGA_", stype, "_20130521_",vtype,"_no_ov_min100");
  return (paste0(name, extension))
}

#load or transform original data
loadOriginalData<-function(stype, vtype){
  name<-getFile(stype, vtype)
  file<-paste0(name, ".t.csv")
  if(!file.exists(file)){
    data<-read.table(paste0(name, ".tsv"), header=T, row.names=1, check.names=F)
    data<-t(data)
    saveData(data, file)
  }else{
    data<-loadData(file)
  }
  return (data)
}

local({
  colLab <<- function(n) {
    if(is.leaf(n)) {
      a <- attributes(n)
      ss <- unlist(strsplit(a$label, "_"))
      col<-cols[ss[1]]
      attr(n, "nodePar") <-
        c(a$nodePar, list(lab.col = as.character(col), lab.font = 1))
    }
    n
  }
})

#function used to draw dendrogram
drawDendrogram<-function(x, filename, width=200000){
  dhc <- as.dendrogram(hc <- hclust(dist(x), "complete"))
  dL <- dendrapply(dhc, colLab)
  
  png(filename, width=width, height=2000,res=300)
  par(mar=c(12, 2, 2, 2))
  plot(dL)
  dev.off()
}

#function used to draw dendrogram
drawHeatmap<-function(x, filename, samplecolor, heatmapcolor, width=200000){
  png(filename, width=width, height=2000,res=300)
  par(mar=c(12, 2, 2, 2))
  heatmap.2(x, ColSideColors = samplecolor, col = heatmapcolor, dendrogram="none" )
  dev.off()
}

#function used to draw headmapplus
drawHeatmapPlus<-function(x, filename, samplecolor, heatmapcolor, width=200000, height=2000, cexRow=1.0, cexCol=1.0, main = NULL, xlab = NULL, ylab = NULL, Colv=NULL){
  if(filename != ""){
    png(filename, width=width, height=height,res=300)
  }
  par(mar=c(12, 2, 2, 2))
  heatmap.plus(x, ColSideColors = samplecolor, col = heatmapcolor,cexRow=cexRow, cexCol=cexCol, main=main, xlab=xlab, ylab=ylab, Colv=Colv)
  if(filename != ""){
    dev.off()
  }
}

#draw by tumor only
drawHeatmapPlusTumor<-function(x, filename, heatmapcolor, width=200000, height=2000, cexRow=1.0, cexCol=1.0, hidesamples=TRUE, hidegenes=FALSE, main = NULL, xlab = NULL, ylab = NULL, Colv = NULL){
  d<-as.matrix(x)
  groupcols<-getSampleColors(d)
  clab  <- matrix(c(rep("white", ncol(d)), groupcols), ncol=2, byrow=FALSE)
  colnames(clab)<-c("", "Tumor")
  if(hidesamples){
    colnames(d)<-c(rep("", ncol(d)))
  }
  if(hidegenes){
    rownames(d)<-c(rep("", nrow(d)))
  }
  drawHeatmapPlus(d, filename, samplecolor=clab, heatmapcolor= heatmapcolor, width, height, cexRow, cexCol,main=main, xlab=xlab, ylab=ylab,Colv=Colv)
}

#draw by tumor and low/high group
drawHeatmapPlusTumorAndExpressionGroup<-function(x, filename, heatmapcolor, width=200000, height=2000, cexRow=1.0, cexCol=1.0, hidesamples=TRUE, hidegenes=FALSE, main = NULL, xlab = NULL, ylab = NULL){
  d<-as.matrix(x)
  groupcols<-getSampleColors(d)
  half<-ncol(d) / 2
  tbcolor<-c(rep("green", half), rep("red", half))
  clab <- matrix(c(tbcolor, groupcols), ncol=2, byrow=FALSE)
  rownames(clab)<-colnames(x)
  colnames(clab)<-c("Low/High","Tumor")
  if(hidesamples){
    colnames(d)<-c(rep("", ncol(d)))
  }
  if(hidegenes){
    rownames(d)<-c(rep("", nrow(d)))
  }
  drawHeatmapPlus(d, filename, samplecolor=clab, heatmapcolor= heatmapcolor, width, height, cexRow, cexCol, main=main, xlab=xlab, ylab=ylab)
}

loadRNAseqData<-function(vtype, xfile){
  if(file.exists(xfile)){
    x<-loadData(xfile)
  }else{
    samples<-read.table("TCGA_mirnaseq_20130521_countexpr_no_ov_min100.t.90.rank.target.lowhigh_no_igg_kirc.sample.tsv", header=T, check.names=F)
    lowfile<-paste0("TCGA_rnaseqv2_20130521_",vtype,"_no_ov_min100_no_igg_kirc_miRNAlow.csv")
    highfile<-paste0("TCGA_rnaseqv2_20130521_",vtype,"_no_ov_min100_no_igg_kirc_miRNAhigh.csv")
    if(!file.exists(lowfile)){
      data<-read.table(paste0("TCGA_rnaseqv2_20130521_",vtype,"_no_ov_min100.tsv"), header=T, row.names=1, check.names=F)
      low<-data[rownames(data) %in% samples$low,]
      low<-t(low)
      saveData(low, lowfile)
      high<-data[rownames(data) %in% samples$high,]
      high<-t(high)
      saveData(high, highfile)
      rm(data)
    }else{
      low<-loadData(lowfile)
      high<-loadData(highfile)
    }
    x<-cbind(low, high)
    saveData(x, xfile)
  }
  return (x)
}

getTargetMirnaRankOrdered<-function(mirna){
  mirnarank<-apply(mirna, 2, function(x) rank(x,ties.method="average"))
  targetMirnaRank<-mirnarank[rownames(mirnarank) %in% targetnames,]
  targetMirnaRankMean<-apply(targetMirnaRank, 2, function(x) mean(x))
  targetMirnaRankOrdered<-rbind(targetMirnaRank, targetMirnaRankMean)
  targetMirnaRankOrdered<-targetMirnaRankOrdered[,order(targetMirnaRankMean)]
  rownames(targetMirnaRankOrdered)[6] <- "meanOfRank"
  return (targetMirnaRankOrdered)
}

#heatmap color
hmcols <- colorRampPalette(c("green", "black", "red"))(256)

targetnames<-c("hsa-mir-141", "hsa-mir-200a", "hsa-mir-200b", "hsa-mir-200c", "hsa-mir-429")