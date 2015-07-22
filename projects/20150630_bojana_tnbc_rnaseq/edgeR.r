
setwd("H:/shengquanhu/projects/Jennifer/20150630_bojana_tnbc/star_deseq2/result")  

showLabelInPCA<-1
showDEGeneCluster<-1
pvalue<-0.05
foldChange<-2
minMedianInGroup<-5

library("edgeR")

data<-read.table("H:/shengquanhu/projects/Jennifer/20150630_bojana_tnbc/star_genetable/result/20150630_bojana_tnbc_gene.count",row.names=1, header=T, check.names=F)

isDataNumeric = unlist(lapply(data[1,], function(x){is.numeric(x)}))
index = 1
while(!all(isDataNumeric[index:ncol(data)])){
  index = index + 1
}

indecies<-c(1:(index-1))
countData<-data[,c(index:ncol(data))]

countData[is.na(countData)] <- 0
countData<-round(countData)

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

#different expression analysis
y<-DGEList(countData, genes=curdata[,1:3])
colnames(y) <- colnames(countData)

condition<-designData$Condition
design <- model.matrix(~condition)
rownames(design)<-colnames(y)

y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)

de <- decideTestsDGE(lrt)
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

cat("EdgeR finished.\n")

select<-(!is.na(res$padj)) & (res$padj<pvalue) & ((res$log2FoldChange >= log2(foldChange)) | (res$log2FoldChange <= -log2(foldChange)))

if(length(indecies) > 0){
  inddata<-curdata[,indecies,drop=F]
  tbb<-cbind(inddata, comparisonData, res)
}else{
  tbb<-cbind(comparisonData, res)
}
tbbselect<-tbb[select,,drop=F]

tbb<-tbb[order(tbb$padj),,drop=F]
write.csv(as.data.frame(tbb),paste0(prefix, "_DESeq2.csv"))

tbbselect<-tbbselect[order(tbbselect$padj),,drop=F]
write.csv(as.data.frame(tbbselect),paste0(prefix, "_DESeq2_sig.csv"))
