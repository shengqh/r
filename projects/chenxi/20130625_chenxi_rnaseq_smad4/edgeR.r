library("edgeR")

glmversion<-0
data<-read.table("H:/shengquanhu/projects/chenxi/20130625_chenxi_rnaseq_smad4/cuffdiff/result/smad4.count", sep="\t", header=T, row.names=1, check.names=F)
data<-data[order(rownames(data)),]
countdata<-data[,3:ncol(data)]
countdata<-countdata[,order(colnames(countdata))]

groups<-read.table("H:/shengquanhu/projects/chenxi/20130625_chenxi_rnaseq_smad4/cuffdiff/result/smad4_group_sample.map", sep="\t", header=T, row.names=1, check.names=F)
groups<-groups[order(groups$SAMPLE_NAME),]

colnames(countdata)<-groups$GROUP_SAMPLE
countdata<-countdata[,order(colnames(countdata))]
groups<-groups[order(groups$GROUP_SAMPLE),]

group<-groups$GROUP

y <- DGEList(counts=countdata,group=group)

if(glmversion){
  design <- model.matrix(~group)
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTrendedDisp(y,design)
  y <- estimateGLMTagwiseDisp(y,design)
  fit <- glmFit(y,design)
  lrt <- glmLRT(fit,coef=2)
  diff<-topTags(lrt, n=nrow(countdata))$table
  diff<-diff[order(rownames(diff)),]
  final<-cbind(data[,1:2], countdata, diff)
  final<-final[order(final$FDR), ]
}else{
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  et<-exactTest(y)
  diff<-topTags(et, n=nrow(countdata))$table
  diff<-diff[order(rownames(diff)),]
  final<-cbind(data[,1:2], countdata, diff)
  final<-final[order(final$FDR), ]
}

write.csv(final, "H:/shengquanhu/projects/chenxi/20130625_chenxi_rnaseq_smad4/cuffdiff/result/smad4.diff.csv")
