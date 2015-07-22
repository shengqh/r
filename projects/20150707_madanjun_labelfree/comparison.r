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
data$padj<-p.adjust(aovpvalues, method="fdr")

sigdata<-data[data$padj < 0.05,]

drawdata<-sigdata[,c(1:(ncol(sigdata)-1))]
  
png(file="ProteinIntensities.sigProteins.heatmap.png", width=4000, height=3000, res=300)
hm<-heatmap3(drawdata, 
             col = hmcols, 
             ColSideColors = gsColors, 
             margins=c(12,5), 
             scale="r", 
             labRow=NA,
             main=paste0("Hierarchical Cluster Using ", nrow(sigdata), " Proteins"),  
             cexCol=cexCol,
             keep.dendro = TRUE,
             legendfun=function() showLegend(legend=as.character(levels(daypoints)), col=tpColors,cex=1.0,x="center"))
dev.off()

png(file="ProteinIntensities.sigProteins.heatmap.png", width=4000, height=3000, res=300)
hm<-heatmap3(drawdata, 
             col = hmcols, 
             ColSideColors = gsColors, 
             margins=c(12,5), 
             scale="r", 
             labRow=NA,
             main=paste0("Hierarchical Cluster Using ", nrow(sigdata), " Proteins"),  
             cexCol=cexCol,
             keep.dendro = TRUE,
             legendfun=function() showLegend(legend=as.character(levels(daypoints)), col=tpColors,cex=1.0,x="center"))
dev.off()

hc.rows<- hclust(as.dist(1 - cor(t(drawdata), use = "pa")))

ct<- cutree(hc.rows, h=1.5)
plot(hc.rows)
rect.hclust(hc.rows, h=1.5)

c1<-ct[ct==4]
c2<-ct[ct==3]
c3<-ct[ct==2]
c4<-ct[ct==1]

ct<-cutree(hc.rows, h=0.75)
plot(hc.rows)
rect.hclust(hc.rows, h=0.75)

c44<-ct[names(ct) %in% names(c4)]
c4<-c44[c44==8]
c5<-c44[c44==2]
c6<-c44[c44==1]

genecluster<-rbind(data.frame(gene=names(c1), cluster="black"),
                   data.frame(gene=names(c2), cluster="red"),
                   data.frame(gene=names(c3), cluster="green"),
                   data.frame(gene=names(c4), cluster="blue"),
                   data.frame(gene=names(c5), cluster="pink"),
                   data.frame(gene=names(c6), cluster="purple")                   )

genecluster$cluster<-as.factor(genecluster$cluster)
rownames(genecluster)<-genecluster$gene

geneColors<-as.matrix(data.frame(Cluster=genecluster[rownames(drawdata),"cluster"]))

png(file="ProteinIntensities.sigProteins.times.heatmap.png", width=4000, height=3000, res=300)
heatmap3(drawdata, 
             Colv=NA,
             col = hmcols, 
             RowSideColors = geneColors, 
             ColSideColors = gsColors, 
             margins=c(12,5), 
             scale="r", 
             labRow=NA,
             main=paste0("Hierarchical Cluster Using ", nrow(sigdata), " Proteins"),  
             cexCol=cexCol,
             keep.dendro = TRUE,
             legendfun=function() showLegend(legend=as.character(levels(daypoints)), col=tpColors,cex=1.0,x="center"))
dev.off()

sigdata$cluster<-genecluster[rownames(sigdata),"cluster"]

write.csv(sigdata,file="ProteinIntensities.sigProteins.cluster.csv")
