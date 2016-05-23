library(reshape2)
library(scales)

setwd("Z:/Shared/Labs/Brown,J/tiger/20160509_brown_3436/star_genetable_deseq2/result")

foldChange<-2
pvalue<-0.05

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

dmso<-read.csv("DMSO_vs_FED_min5_DESeq2.csv", header=T, stringsAsFactors=F)
dmso<-dmso[,c("Feature_gene_name", "baseMean", "log2FoldChange", "padj" )]
dmso<-dmso[which(dmso$padj<=pvalue & abs(dmso$log2FoldChange)>=log2(foldChange)),]
dmso$log10BaseMean<-log10(dmso$baseMean)
dmso$Sample<-"DMSO"

jq1<-read.csv("JQ1_vs_FED_min5_DESeq2.csv", header=T, stringsAsFactors=F)
jq1<-jq1[,c("Feature_gene_name", "baseMean", "log2FoldChange", "padj" )]
jq1<-jq1[which(jq1$padj<=pvalue & abs(jq1$log2FoldChange)>=log2(foldChange)),]
jq1$log10BaseMean<-log10(jq1$baseMean)
jq1$Sample<-"JQ1"

jq1<-jq1[jq1$Feature_gene_name %in% dmso$Feature_gene_name,]
dmso<-dmso[dmso$Feature_gene_name %in% jq1$Feature_gene_name,]

changeColours<-c(DMSO="blue",JQ1="red")
diffResult<-rbind(dmso, jq1)

png(filename="DMSO_JQ1_volcanoPlot.png", width=3000, height=3000, res=300)
p<-ggplot(diffResult,aes(x=log2FoldChange,y=padj))+
  geom_point(aes(colour=Sample), size=2)+
  scale_color_manual(values=changeColours,guide = FALSE)+
  scale_y_continuous(trans=reverselog_trans(10),name=bquote(Adjusted~p~value))+
  scale_x_continuous(name=bquote(log[2]~Fold~Change))+
  geom_hline(yintercept = 1,colour="grey",linetype = "dotted")+
  geom_vline(xintercept = 0,colour="grey",linetype = "dotted")+
  theme_bw()+
  theme(axis.text = element_text(colour = "black",size=30),
        axis.title = element_text(size=30))
print(p)
dev.off()
