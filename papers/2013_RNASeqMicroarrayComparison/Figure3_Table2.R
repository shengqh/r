library("grid")
library("VennDiagram")

rootdir<-"E:/sqh/Dropbox/papers/2013-RNAseqMicroarrayComparison";

tiff(filename=paste0(rootdir,"/Figure3_RelativeComparison.tif"), width=3000, height=1500, res=300, compression="lzw")

par(mfrow=c(1,2))

fx=0.07
fy=0.96
fy2=1-(1-fy)*2
fsize=16

#relative comparison
par(mar=c(14,9,3.2,3.5))

datafile<-paste0(rootdir, "/Figure3_BrcaRelativeComparison.tsv");
data<-read.table(datafile, header=T,row.names=1, check.names=F)
labels=rep("",length(colnames(data)))
bxp <- boxplot(data, names=labels, ylab="Spearman correlation\ncoefficient", ylim=c(0.25, 1.0), medcol="red", boxwex=0.45)

text(x=c(1:3) + 0.42, y = bxp$stats[3,1:3], labels=bxp$stats[3,1:3], cex=0.5)
text(x=c(1:3), y=rep(0.1,3), labels = colnames(data), xpd = TRUE,cex=0.7)

lines(x=c(1,1,2,2),y=c(0.90,0.92,0.92,0.90))
text(x=1.5, y=0.95, labels=c('** p < 0.001'), cex=0.5)

wilcox.test(data[,1], data[,2], paired=TRUE)

#venndiagram
datafile<-paste0(rootdir, "/Figure3_BRCA_PValue.csv");
data<-read.csv(datafile, header=T, row.names=1, check.names=F)

agi<-subset(data, (abs(data[,1]) >= 1) & (data[,3] <= 0.01))
write.csv(agi, paste0(rootdir, "/Figure3_Agilent_significant.csv"))

rnaseq<-subset(data, (!is.na(data[,6])) & (abs(data[,4]) >= 1) & (data[,6] <= 0.01))
write.csv(rnaseq, paste0(rootdir, "/Figure3_RNAseq_significant.csv"))

common<-subset(data, (abs(data[,1]) >= 1) & (data[,3] <= 0.01) & (!is.na(data[,6])) & (abs(data[,4]) >= 1) & (data[,6] <= 0.01))
union<-subset(data, ((abs(data[,1]) >= 1) & (data[,3] <= 0.01)) || ((!is.na(data[,6])) & (abs(data[,4]) >= 1) & (data[,6] <= 0.01)))

pushViewport(viewport(x=0, width = 0.5, y=0.75, height=0.5, just="left"))
grid.text("a",x=fx, y=fy2, gp=gpar(fontsize=fsize))
upViewport()
pushViewport(viewport(x=0, width = 0.5, y=0.25, height=0.5, just="left"))
glist<-draw.pairwise.venn(
  area1 = nrow(rnaseq),
  area2 = nrow(agi),
  cross.area = nrow(common),
  category = c("", ""),
  fill = c("blue", "red"),
  lty = "blank",
  lwd = 0.3,
  cex = 1,
  cat.cex = 1,
  cat.pos = c(285, 108),
  cat.dist = c(0.11, 0.09),
  cat.just = list(c(-1, -1), c(1, 1)),
  ext.pos = 30,
  ext.dist = -0.05,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.line.lty = "dashed",
  ind=F)
tmp <- 
glist <- gList(glist, 
               textGrob(label = "RNAseq", x = 0.18, y = 0.78, just = c(0,0), gp = gpar(col = "blue", cex = 1, fontface = "plain", fontfamily = "serif")),
               textGrob(label = "Agilent", x = 0.71, y = 0.78, just = c(0,0), gp = gpar(col = "red", cex = 1, fontface = "plain", fontfamily = "serif")))
grid.draw(glist)
grid.text("b",x=fx, y=fy2, gp=gpar(fontsize=fsize))
upViewport()

pushViewport(viewport(x = 0.5, width = 0.5, just="left"))
grid.text("c",x=0.1, y=fy, gp=gpar(fontsize=fsize))
upViewport()


#significant overlap
par(mar=c(6,6,3,1))

cols<-apply(data, 1, function(x) {
  agisig<-(abs(x["AgilentFold"]) >= 1) & (x["AgilentFDR"] <= 0.01)
  rpkmsig<-!is.na(x["RPKMFDR"]) & (abs(x["RPKMFold"]) >= 1) & (x["RPKMFDR"] <= 0.01)
  
  if(agisig & rpkmsig){
    if( ((x["AgilentFold"] < 0) & (x["RPKMFold"] < 0) ) || ((x["AgilentFold"] > 0) & (x["RPKMFold"] > 0))){
      return("green2")
    }else{
      return("red")
    }
  }
  else if(agisig){
    return ("yellow3")
  }
  else if(rpkmsig){
    return ("blue")
  }
  else{
    return ("black")
  }
})

coltabl<-table(cols)
inconsistent<-coltabl["red"] * 100 / (coltabl["red"] + coltabl["green2"])

inconsistentData<-subset(data, cols=="red")[,c(1:6)]
colnames(inconsistentData)<-c("log2(AgilentFoldChange)",  "pValue(Agilent)",	"fdr(Agilent)",	"log2(RPKMFoldChange)",	"pValue(RPKM)",	"fdr(RPKM)")
write.csv(inconsistentData, paste0(rootdir, "/Table2_DiffOritation.csv"))

plot(data$RPKMFold, data$AgilentFold,  ylab="log2(Agilent fold change)", xlab="log2(RNAseq fold change)", ylim=c(-10,10), col=cols, cex=0.2)

legcols<-c("red","green2", "blue", "yellow3", "black")
legend("topleft", c(sprintf("Inconsistent (%0.1f%%)", inconsistent),
                    sprintf("Consistent (%0.1f%%)", 100 - inconsistent),
                    "RNAseq only", "Microarray only", "Neither"), pch=c(1,1), col=legcols, bty="n", text.col= legcols)

dev.off();

datafile<-paste0(rootdir, "/Figure3_BrcaRelativeComparison.tsv");
data<-read.table(datafile, header=T,row.names=1, check.names=F)
aa<-apply(data,2, function(x) c(median(x), min(x), max(x)))
rownames(aa)<-c("median","minimum","maximum")
write.csv(aa, paste0(rootdir, "/Figure3_BrcaRelativeComparison.summary.csv"))