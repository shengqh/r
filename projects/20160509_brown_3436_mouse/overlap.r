library(reshape2)

setwd("Z:/Shared/Labs/Brown,J/tiger/20160509_brown_3436/star_genetable_deseq2/result")

comparisons = c("DMSO_vs_FED", "JQ1_vs_FED", "CAPTISOL_vs_FED", "dBET_vs_FED", "JQ1_vs_DMSO", "dBET_vs_CAPTISOL")

res=NULL
comp=comparisons[2]
for (comp in comparisons){
  compfile=paste0(comp, "_min5_DESeq2.csv")
  compdata<-read.csv(compfile, header=T, stringsAsFactors=F)
  compsigdata<-compdata[compdata$padj<=0.05 & abs(compdata$log2FoldChange)>=1,]
  compsigdata$Comparison=comp
  compsigdata<-compsigdata[,c("Comparison", "Feature_gene_name", "FoldChange" )]
  if(is.null(res)){
    res=compsigdata
  }else{
    res=rbind(res,compsigdata)
  }
}

dres=dcast(res, Feature_gene_name ~ Comparison,fun.aggregate=median)
dres[is.na((dres))]=""
write.csv(dres, file="sig_overlap.csv", row.names=F)
