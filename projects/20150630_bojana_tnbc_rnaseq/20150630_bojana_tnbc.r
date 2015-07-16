
setwd("H:/shengquanhu/projects/Jennifer/20150630_bojana_tnbc/star_deseq2/result")  

data<-read.table("H:/shengquanhu/projects/Jennifer/20150630_bojana_tnbc/star_genetable/result/20150630_bojana_tnbc_gene.count",row.names=1, header=T, check.names=F)

showLabelInPCA<-1
showDEGeneCluster<-1
pvalue<-0.05
foldChange<-2
minMedianInGroup<-5

library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")

##Solving node stack overflow problem start###
#when there are too many genes, drawing dendrogram may failed due to node stack overflow,
#It could be solved by forcing stats:::plotNode to be run as interpreted code rather then byte-compiled code via a nasty hack.
#http://stackoverflow.com/questions/16559250/error-in-heatmap-2-gplots/25877485#25877485

# Convert a byte-compiled function to an interpreted-code function 
unByteCode <- function(fun)
{
  FUN <- eval(parse(text=deparse(fun)))
  environment(FUN) <- environment(fun)
  FUN
}

# Replace function definition inside of a locked environment **HACK** 
assignEdgewise <- function(name, env, value)
{
  unlockBinding(name, env=env)
  assign( name, envir=env, value=value)
  lockBinding(name, env=env)
  invisible(value)
}

# Replace byte-compiled function in a locked environment with an interpreted-code
# function
unByteCodeAssign <- function(fun)
{
  name <- gsub('^.*::+','', deparse(substitute(fun)))
  FUN <- unByteCode(fun)
  retval <- assignEdgewise(name=name,
                           env=environment(FUN),
                           value=FUN
  )
  invisible(retval)
}

# Use the above functions to convert stats:::plotNode to interpreted-code:
unByteCodeAssign(stats:::plotNode)

# Now raise the interpreted code recursion limit (you may need to adjust this,
#  decreasing if it uses to much memory, increasing if you get a recursion depth error ).
options(expressions=5e4)

##Solving node stack overflow problem end###

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

drawHCA<-function(prefix, rldselect, designData, conditionColors, gnames){
  htfile<-paste0(prefix, "_DESeq2-vsd-heatmap.png")
  cat("saving HCA to ", htfile, "\n")
  genecount<-nrow(rldselect)
  if(genecount > 2){
    png(filename=htfile, width=3000, height =3000, res=300)
    cexCol = max(1.0, 0.2 + 1/log10(ncol(rldselect)))
    
    subtypeColors<-rainbow(length(levels(designData$Subtype)))
    treatmentColors<-rainbow(length(levels(designData$Treatment)))
    
    gsColors<-as.matrix(data.frame(Response=conditionColors, 
                                   Subtype=subtypeColors[designData$Subtype], 
                                   Treatment=treatmentColors[designData$Treatment]))
    
    legendnames<-c(paste("Response:", as.character(levels(designData$Response))), 
                   "",
                   paste("Subtype:", as.character(levels(designData$Subtype))), 
                   "",
                   paste("Treatment:",as.character(levels(designData$Treatment))))
    legendcolors<-c("green", "blue", "red", "white", subtypeColors, "white", treatmentColors )
    
    heatmap3(rldselect, 
             col = hmcols, 
             ColSideColors = gsColors, 
             margins=c(12,5), 
             scale="r", 
             dist=dist, 
             labRow=NA,
             main=paste0("Hierarchical Cluster Using ", genecount, " Genes"),  
             cexCol=cexCol, 
             legendfun=function() showLegend(legend=legendnames, col=legendcolors,cex=1.0,x="center"))
    dev.off()
  }
}

drawPCA<-function(prefix, rldmatrix, showLabelInPCA, designData, conditionColors){
  filename<-paste0(prefix, "_DESeq2-vsd-pca.png")
  cat("saving PCA to ", filename, "\n")
  png(filename=filename, width=3000, height=3000, res=300)
  pca<-prcomp(t(rldmatrix))
  supca<-summary(pca)$importance
  pcadata<-data.frame(pca$x)
  pcalabs=paste0(colnames(pcadata), "(", round(supca[2,] * 100), "%)")
  pcadata["sample"]<-row.names(pcadata)
  
  if(showLabelInPCA){
    g <- ggplot(pcadata, aes(x=PC1, y=PC2, label=sample)) + 
      geom_text(vjust=-0.6, size=4) +
      geom_point(col=conditionColors, size=4) + 
      scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) +
      scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
      geom_hline(aes(0), size=.2) + 
      geom_vline(aes(0), size=.2) + 
      xlab(pcalabs[1]) + ylab(pcalabs[2])
  }else{
    g <- ggplot(pcadata, aes(x=PC1, y=PC2)) + 
      geom_point(col=conditionColors, size=4) + 
      labs(color = "Group") +
      scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) + 
      scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
      geom_hline(aes(0), size=.2) + 
      geom_vline(aes(0), size=.2) +
      xlab(pcalabs[1]) + ylab(pcalabs[2]) + 
      theme(legend.position="top")
  }
  
  print(g)
  dev.off()
}

isDataNumeric = unlist(lapply(data[1,], function(x){is.numeric(x)}))
index = 1
while(!all(isDataNumeric[index:ncol(data)])){
  index = index + 1
}

indecies<-c(1:(index-1))

designData<-read.csv("design.csv")
countData<-data[,as.character(designData$Sample)]
countData[is.na(countData)] <- 0
countData<-round(countData)

designData$Condition<-as.factor(ifelse(designData$Response=="Group 3","Group 3","Group 1+2"))

comparisonName<-"ClinicalResponse"
prefix<-comparisonName
curdata<-data
if(minMedianInGroup > 0){
  conds<-as.character(levels(designData$Response))
  data1<-countData[, colnames(countData) %in% designData$Sample[designData$Response==conds[1]]]
  data2<-countData[, colnames(countData) %in% designData$Sample[designData$Response==conds[2]]]
  data3<-countData[, colnames(countData) %in% designData$Sample[designData$Response==conds[3]]]
  med1<-apply(data1, 1, median) >= minMedianInGroup
  med2<-apply(data2, 1, median) >= minMedianInGroup
  med3<-apply(data3, 1, median) >= minMedianInGroup
  
  med<-med1 | med2 | med3

  countData<-countData[med,]
  cat(nrow(countData), " genes with minimum median count in group larger or equals than ", minMedianInGroup, "\n")
  prefix<-paste0(comparisonName, "_min", minMedianInGroup)
  curdata<-data[med,]
}

notEmptyData<-apply(countData, 1, max) > 0
countData<-countData[notEmptyData,]
curdata<-curdata[notEmptyData,]

rownames(designData)<-designData$Sample
conditionColors<-as.matrix(data.frame(Response=c("green", "blue", "red")[designData$Response]))

#some basic graph
dds=DESeqDataSetFromMatrix(countData = countData,
                           colData = designData,
                           design = ~1)

colnames(dds)<-colnames(countData)

#draw density graph
rldmatrix<-as.matrix(log2(counts(dds,normalized=FALSE) + 1))
rsdata<-melt(rldmatrix)
colnames(rsdata)<-c("Gene", "Sample", "log2Count")
png(filename=paste0(prefix, "_DESeq2-log2-density.png"), width=6000, height=4000, res=300)
g<-ggplot(rsdata) + geom_density(aes(x=log2Count, colour=Sample)) + xlab("DESeq2 log2 transformed count")
print(g)
dev.off()

png(filename=paste0(prefix, "_DESeq2-log2-density-individual.png"), width=6000, height=4000, res=300)
g<-ggplot(rsdata) + geom_density(aes(x=log2Count, colour=Sample)) + facet_wrap(~Sample, scales = "free") + xlab("DESeq2 log2 transformed count")
print(g)
dev.off()

#varianceStabilizingTransformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
assayvsd<-assay(vsd)
write.csv(assayvsd, file=paste0(prefix, "_DESeq2-vsd.csv"))

vsdiqr<-apply(assayvsd, 1, IQR)
assayvsd<-assayvsd[order(vsdiqr, decreasing=T),]

rldmatrix=as.matrix(assayvsd)

#draw pca graph
drawPCA(paste0(prefix,"_geneAll"), rldmatrix, showLabelInPCA, designData, conditionColors)

#draw heatmap
drawHCA(paste0(prefix,"_gene500"), rldmatrix[1:min(500, nrow(rldmatrix)),,drop=F], designData, conditionColors)

drawHCA(paste0(prefix,"_geneAll"), rldmatrix, designData, conditionColors)

dds=DESeqDataSetFromMatrix(countData = countData,
                           colData = designData,
                           design = ~ Condition)

dds <- DESeq(dds)
res<-results(dds,cooksCutoff=FALSE)

cat("DESeq2 finished.\n")

select<-(!is.na(res$padj)) & (res$padj<pvalue) & ((res$log2FoldChange >= log2(foldChange)) | (res$log2FoldChange <= -log2(foldChange)))

if(length(indecies) > 0){
  inddata<-curdata[,indecies,drop=F]
  tbb<-cbind(inddata, countData, res)
}else{
  tbb<-cbind(countData, res)
}
tbbselect<-tbb[select,,drop=F]

tbb<-tbb[order(tbb$padj),,drop=F]
write.csv(as.data.frame(tbb),paste0(prefix, "_DESeq2.csv"))

tbbselect<-tbbselect[order(tbbselect$padj),,drop=F]
write.csv(as.data.frame(tbbselect),paste0(prefix, "_DESeq2_sig.csv"))

if(showDEGeneCluster){
  siggenes<-rownames(rldmatrix) %in% rownames(tbbselect)
  
  #nonDEmatrix<-rldmatrix[!siggenes,,drop=F]
  #drawPCA(paste0(prefix,"_geneNotDE"), nonDEmatrix, showLabelInPCA, designData, conditionColors)
  #drawHCA(paste0(prefix,"_geneNotDE"), nonDEmatrix, designData, conditionColors)
  
  DEmatrix<-rldmatrix[siggenes,,drop=F]
  drawHCA(paste0(prefix,"_geneDE"),DEmatrix, designData, conditionColors)
}

