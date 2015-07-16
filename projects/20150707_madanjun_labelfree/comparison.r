library("heatmap3")

setwd("E:/shengquanhu/projects/20150707_madanjun_labelfree")

rawdata<-read.table("ProteinIntensities.tsv", header=T, sep="\t", row.names=1)
data<-log2(rawdata[,c(2:17)])

daypoints<-paste("Dayl", gsub(".*Low_(\\d+).*\\d", "\\1", colnames(data)))
reppoints<-paste("rep", gsub(".*(\\d)$", "\\1", colnames(data)))
colnames(data)<-paste(daypoints, reppoints)

daypoints<-factor(daypoints, levels=unique(daypoints))
hmcols <- colorRampPalette(c("green", "black", "red"))(256)
tpColors<-rainbow(length(unique(daypoints)))
gsColors<-as.matrix(data.frame(Day=tpColors[daypoints]))
cexCol = max(1.0, 0.2 + 1/log10(ncol(data)))

png(file="ProteinIntensities.allProteins.heatmap.png", width=4000, height=3000, res=300)
heatmap3(data, 
         col = hmcols, 
         ColSideColors = gsColors, 
         margins=c(12,5), 
         scale="r", 
         dist=dist, 
         labRow=NA,
         main=paste0("Hierarchical Cluster Using ", nrow(data), " Proteins"),  
         cexCol=cexCol,
         legendfun=function() showLegend(legend=as.character(levels(daypoints)), col=tpColors,cex=1.0,x="center"))
dev.off()

aovpvalues<-apply(data,1,function(x){
  df<-data.frame(Intensity=unlist(x), Group=daypoints)
  fit <- aov(Intensity ~ Group, data=df)
  summary(fit)[[1]][["Pr(>F)"]][[1]]
})
aovpadj<-p.adjust(aovpvalues, method="fdr")
sigdata<-data[aovpadj < 0.05,]

write.csv(sigdata, file="anova_sig.csv")

png(file="ProteinIntensities.sigProteins.heatmap.png", width=4000, height=3000, res=300)
heatmap3(sigdata, 
         col = hmcols, 
         ColSideColors = gsColors, 
         margins=c(12,5), 
         scale="r", 
         dist=dist, 
         labRow=NA,
         main=paste0("Hierarchical Cluster Using ", nrow(sigdata), " Proteins"),  
         cexCol=cexCol,
         legendfun=function() showLegend(legend=as.character(levels(daypoints)), col=tpColors,cex=1.0,x="center"))
dev.off()