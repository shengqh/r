setwd("H:/shengquanhu/projects/Jennifer/20150805_johanna_neuroblastoma")
data<-read.table("CCLE_CellLines.tsv", sep="\t", header=T)

defile<-"CCLE_CellLines_DE.tsv"
if(!file.exists(defile)){
  countdata<-data[,-1]
  
  library(limma)
  
  combn<-factor(c(rep("MYCN_HIGH",11), rep("MYCN_LOW",4)))
  design <- model.matrix(~combn) 
  
  fit <- lmFit(countdata, design)
  efit <- eBayes(fit)
  tp<-topTable(efit, coef=2, number=nrow(countdata))
  tp<-tp[tp$adj.P.Val < 0.05,]
  tp<-tp[abs(tp$logFC) >= 1,]
  
  tpnames<-rownames(tp)
  
  finaldata<-cbind(data.frame(Gene=data[tpnames,"Description"]), tp[,c(1,4,5)], countdata[tpnames,])
  write.table(finaldata, file=defile, row.names=F, sep="\t")
}

expectdata<-data[data$Description %in% c("MYC","MYCN", "MYCL1","TFPI2"),]
meltdata<-melt(expectdata)
colnames(meltdata)<-c("Gene","Sample","RMA")

png(file="CCLE_MYC_MYCN_MYCL1_TFPI2.png", height=1300, width=max(2000, 50 * (ncol(expectdata) - 1)), res=300)
g<-ggplot(meltdata, aes(x=Sample, y=RMA))+ 
  geom_point(aes(colour=Gene), size=3) +
  ggtitle("CCLE MYC/MYCN/MYCL1/TFPI2 Expression Values") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10, hjust=0, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        plot.title=element_text(face="bold", size=16))
print(g)
dev.off()