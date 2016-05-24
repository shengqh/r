library(reshape2)
library(scales)
library(ggplot2)
library(cowplot)

setwd("E:/projects/20160509_brown_3436/star_genetable_deseq2/result")

foldChange<-2
pvalue<-0.05

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

dmso<-read.csv("DMSO_vs_FED_min5_DESeq2_sig.csv", header=T, stringsAsFactors=F, row.names=1)
dmso<-dmso[,c("Feature_gene_name", "baseMean", "log2FoldChange", "padj" )]
dmso<-dmso[which(dmso$padj<=pvalue & abs(dmso$log2FoldChange)>=log2(foldChange)),]
dmso$log10BaseMean<-log10(dmso$baseMean)
dmso$Sample<-"DMSO"

jq1<-read.csv("JQ1_vs_FED_min5_DESeq2_sig.csv", header=T, stringsAsFactors=F, row.names=1)
jq1<-jq1[,c("Feature_gene_name", "baseMean", "log2FoldChange", "padj" )]
jq1<-jq1[which(jq1$padj<=pvalue & abs(jq1$log2FoldChange)>=log2(foldChange)),]
jq1$log10BaseMean<-log10(jq1$baseMean)
jq1$Sample<-"JQ1"

commonNames<-rownames(jq1)[rownames(jq1) %in% rownames(dmso)]
jq1<-jq1[commonNames,]
dmso<-dmso[commonNames,]

sameDirection<-unlist(lapply(commonNames, function(x){
  if(jq1[x, "log2FoldChange"] > 0 & dmso[x, "log2FoldChange"] > 0)
    return ("SAME")
  if(jq1[x, "log2FoldChange"] < 0 & dmso[x, "log2FoldChange"] < 0)
    return ("SAME")
  return ("OPPOSITE")
}))

jq1$Direction<-sameDirection
dmso$Direction<-sameDirection

changeColours<-c("SAME"="blue","OPPOSITE"="red")
diffResult<-rbind(dmso, jq1)

png(filename="DMSO_JQ1_volcanoPlot.png", width=3000, height=3000, res=300)
p<-ggplot(diffResult,aes(x=log2FoldChange,y=padj))+
  geom_point(data=diffResult[diffResult$Direction=="SAME",], aes(colour=Direction), size=2)+
  geom_point(data=diffResult[diffResult$Direction=="OPPOSITE",],aes(colour=Direction), size=2)+
  scale_color_manual(values=changeColours)+
  scale_y_continuous(trans=reverselog_trans(10),name=bquote(Adjusted~p~value))+
  scale_x_continuous(name=bquote(log[2]~Fold~Change))+
  geom_hline(yintercept = 1,colour="grey",linetype = "dotted")+
  geom_vline(xintercept = 0,colour="grey",linetype = "dotted")+
  facet_grid(~Sample)+
  theme_cowplot()+
  theme(axis.text = element_text(colour = "black",size=20),
        axis.title = element_text(size=20),
        legend.text = element_text(colour = "black",size=15),
        legend.title = element_text(size=15))
  print(p)
dev.off()
