library("heatmap3")

setwd("E:/shengquanhu/projects/20150707_madanjun_labelfree")

rawdata<-read.table("ProteinIntensities.tsv", header=T, sep="\t", row.names=1)
data<-log2(rawdata[,c(2:17)])

daypoints<-paste("Dayl", gsub(".*Low_(\\d+).*\\d", "\\1", colnames(data)))
reppoints<-paste("rep", gsub(".*(\\d)$", "\\1", colnames(data)))
colnames(data)<-paste(daypoints, reppoints)

samples<-colnames(data)

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

data$maxLog2FoldChange<-apply(data,1,function(x){
  ad<-2 ** unlist(x)
  fc1<-abs(log2(mean(ad[1:4]) / mean(ad[5:8])))
  fc2<-abs(log2(mean(ad[1:4]) / mean(ad[9:12])))
  fc3<-abs(log2(mean(ad[1:4]) / mean(ad[13:16])))
  fc4<-abs(log2(mean(ad[5:8]) / mean(ad[9:12])))
  fc5<-abs(log2(mean(ad[5:8]) / mean(ad[13:16])))
  fc6<-abs(log2(mean(ad[9:12]) / mean(ad[13:16])))
  max(fc1,fc2,fc3,fc4,fc5,fc6)
})

write.csv(data,file="ProteinIntensities.allProteins.csv")

sigdata<-data[(data$padj < 0.05),]
drawdata<-sigdata[,samples]

png(file="ProteinIntensities.sigProteins.heatmap.png", width=4000, height=3000, res=300)
hm<-heatmap3(drawdata, 
             col = hmcols,
             Colv=c(1:16),
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

png(file="ProteinIntensities.sigProteins.heatmap.cluster.png", width=4000, height=3000, res=300)
heatmap3(drawdata, 
         Colv=c(1:16),
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

sigdata<-sigdata[sigdata$maxLog2FoldChange >= 1,]
drawdata<-sigdata[,samples]

png(file="ProteinIntensities.sigProteins.foldchange2.heatmap.png", width=4000, height=3000, res=300)
heatmap3(drawdata, 
         col = hmcols, 
         Colv = c(1:16),
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

ct<- cutree(hc.rows, h=1.7)
plot(hc.rows)
rect.hclust(hc.rows, h=1.7)

c1<-ct[ct==3]
c2<-ct[ct==2]
c3<-ct[ct==1]

genecluster<-rbind(data.frame(gene=names(c1), cluster="black"),
                   data.frame(gene=names(c2), cluster="pink"),
                   data.frame(gene=names(c3), cluster="Green"))

genecluster$cluster<-as.factor(genecluster$cluster)
rownames(genecluster)<-genecluster$gene

geneColors<-as.matrix(data.frame(Cluster=genecluster[rownames(drawdata),"cluster"]))

png(file="ProteinIntensities.sigProteins.foldchange2.heatmap.cluster.png", width=4000, height=3000, res=300)
heatmap3(drawdata, 
             col = hmcols, 
             Colv = c(1:16),
             RowSideColors=geneColors,
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
write.csv(sigdata,file="ProteinIntensities.sigProteins.foldchange2.cluster.csv")

# 
# stressResponse<-c(0,0,0,0,3,3,3,3,2,2,2,2,1,1,1,1)
# initiateRegeneration<-c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
# dedifferentiation<-c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1)
# unknown1<-c(2,2,2,2,2,2,2,2,1,1,1,1,0,0,0,0)
# unknown2<-c(1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0)
# 
# biocluster<-apply(drawdata,1,function(x){
#   p1<-cor.test(x,stressResponse)$p.value
#   p2<-cor.test(x,initiateRegeneration)$p.value
#   p3<-cor.test(x,dedifferentiation)$p.value
#   p4<-cor.test(x,unknown1)$p.value
#   p5<-cor.test(x,unknown2)$p.value
#   pmin<-min(p1,p2,p3,p4,p5)
#   if(pmin == p1){
#     c("stress response", pmin)
#   }else if(pmin == p2){
#     c("initiate regeneration", pmin)
#   }else if(pmin == p3){
#     c("dedifferentiation", pmin)
#   }else if(pmin == p4){
#     c("unknown1", pmin)
#   }else {
#     c("unknown2", pmin)
#   }
# })
# 
# padjTemplate<-p.adjust(biocluster[2,], method="fdr")
# sigdata$padjTemplate<-padjTemplate
# sigdata$Template<-biocluster[1,]
# write.csv(sigdata,file="ProteinIntensities.sigProteins.foldchange2.biocluster.csv")
