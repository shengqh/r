#Direct comparison between Microarray/RPKM/RSEM
library(grid)
library(ggplot2)
library(reshape)
library(gridExtra)

rootdir<-"E:/sqh/Dropbox/papers/2013-RNAseqMicroarrayComparison";

data<-read.table(paste0(rootdir,"/Figure2_Microarray_RNASeq_Comparison.tsv"), header=T,row.names=1, check.names=F)

tiff(file=paste0(rootdir,"/Figure2_MicroarrayRNASeqComparison.tif"), width=2000, height=2000, res=300, compression="lzw")

grid.newpage()

textx<-element_text(angle=45, hjust=1, color="black", size=6)
titley<-element_text(color="black", size=8)

rnaseq <- ggplot(melt(data[,29:ncol(data)])) + geom_boxplot(aes(x = variable, y = value)) +
  ylim(0.8, 1.0) + 
  ylab("Spearman correlation coefficient") + 
  theme(axis.title.y=titley,
        axis.title.x=element_blank(),
        axis.text.x=textx, 
        axis.ticks = element_blank(),
        plot.margin=unit(c(0.2,0.6,0.2,0.4), "cm"))

affy_rnaseq <- ggplot(melt(data[,24:28])) + geom_boxplot(aes(x = variable, y = value)) +
  ylim(0.5, 1.0) +  
  ylab("Spearman correlation coefficient") + 
  theme(axis.title.y=titley,
        axis.title.x=element_blank(),
        axis.text.x=textx, 
        axis.ticks = element_blank(),
        plot.margin=unit(c(0.2,0.2,0.2,0.6), "cm"))

agilent_rnaseq <- ggplot(melt(data[,4:23])) + geom_boxplot(aes(x = variable, y = value)) + 
  ylim(-0.4, 0.6) + 
  ylab("Spearman correlation coefficient") + 
  theme(axis.title.y=titley,
        axis.title.x=element_blank(),
        axis.text.x=textx, 
        axis.ticks = element_blank(),
        plot.margin=unit(c(0.5,0.6,0.2,0.6), "cm"))

microarray <- ggplot(melt(data[,1:3])) + geom_boxplot(aes(x = variable, y = value)) + 
  ylim(-0.4, 0.6) + 
  ylab("Spearman correlation coefficient") + 
  theme(axis.title.y=titley,
        axis.title.x=element_blank(),
        axis.text.x=textx, 
        axis.ticks = element_blank(),
        plot.margin=unit(c(0.5,0.6,0.2,0.6), "cm"))


pairedfilter <- function(x) {
  if(is.na(x[1]) && is.na(x[2])){
    return (c(NA,NA))
  }

  if(is.na(x[3]) && is.na(x[4]) && is.na(x[5])){
    return (c(NA,NA))
  }
  
  return (c(max(x[1:2], na.rm=TRUE),max(x[3:5], na.rm=TRUE)))
}


grid.arrange(arrangeGrob(rnaseq, affy_rnaseq, ncol=2, widths=c(2,1)), arrangeGrob(agilent_rnaseq, microarray, ncol=2, widths=c(3,1)), nrow=2)

glist <- gList(textGrob(label = "a", x = 0.02, y = 0.95, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif")),
               textGrob(label = "b", x = 0.69, y = 0.95, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif")),
               textGrob(label = "c", x = 0.02, y = 0.47, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif")),
               textGrob(label = "d", x = 0.76, y = 0.47, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif")))
grid.draw(glist)

xadd=0.66
grid.lines(x=c(0.134, 0.179) + xadd, y = c(0.91, 0.91))
grid.lines(x=c(0.220, 0.304) + xadd, y = c(0.91, 0.91))
grid.lines(x=c(0.156, 0.156, 0.262, 0.262) + xadd, y = c(0.91, 0.928, 0.928,0.91))
grid.text("** p < 0.001", x=unit(0.213 + xadd, "npc"), y = unit(0.94, "npc"), gp=gpar(cex=0.6))

dev.off()

affymatrix=apply(data[,24:28],1,function(x) pairedfilter(x))
affymatrix=affymatrix[, !is.na(affymatrix[1,])]
wt<-wilcox.test(x=affymatrix[1,], y=affymatrix[2,], paired=TRUE)
print (paste0("commonsample=",ncol(affymatrix), ", wilcox.pvalue=",wt$p.value))

bxp<-boxplot(data, plot=FALSE)
stats<-rbind(c(min(bxp$stats[3,1:3]), mean(bxp$stats[3,1:3]), max(bxp$stats[3,1:3])),
             c( min(bxp$stats[3,4:12]), mean(bxp$stats[3,4:12]), max(bxp$stats[3,4:12])),
             c( min(bxp$stats[3,13:23]), mean(bxp$stats[3,13:23]), max(bxp$stats[3,13:23])),
             c( min(bxp$stats[3,24:25]), mean(bxp$stats[3,24:25]), max(bxp$stats[3,24:25])),
             c( min(bxp$stats[3,26:28]), mean(bxp$stats[3,26:28]), max(bxp$stats[3,26:28])),
             c( min(bxp$stats[3,29:ncol(data)]), mean(bxp$stats[3,29:ncol(data)]), max(bxp$stats[3,29:ncol(data)])))
rownames(stats)<-c("agilent/affy","agilent/RPKM","agilent/RSEM","affy/RPKM","affy/RSEM","RPKM/RSEM")
colnames(stats)<-c("min(aveS)","mean(aveS)","max(aveS)")
write.csv(stats, file=paste0(rootdir,"/Figure2_MicroarrayRNASeqComparison.stat.csv"))