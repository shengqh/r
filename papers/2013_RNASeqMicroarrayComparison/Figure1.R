library(grid)

myhist <- function(a, xlab, maxvalue, breakcount){
  step = maxvalue / breakcount
  breaks=c(step * (0:breakcount))
  a[a>maxvalue]<-maxvalue;
  hist(a, xlab=xlab, breaks=breaks, main="")
}

printrange<-function(a, title){
  b<-a[!is.na(a)]
  zeros<-b[b==0]
  percentages<-length(zeros) * 100.0/length(b)
  cat(paste("range of", title, "=[", min(b), "-", max(b), "], zeros=", percentages ));
}
 
rootdir<-"E:/sqh/Dropbox/papers/2013-RNAseqMicroarrayComparison";
breakcount=40;
 
datafile<-paste0(rootdir, "/Figure1_Distribution.tsv");
data<-read.table(datafile, header=T, row.names=1, check.names=F, sep="\t")
tiff(file=paste0(rootdir,"/Figure1_Distribution.tif"), width=2000, height=2000, res=300, compression="lzw")
par(mar=c(5,5,2,1))
par(mfrow=c(2,2))
xper=-0.3
yper=1.1
hist(data[,1], xlab="Affymetrix RMA value",breaks=c(0.4*(0:40)), main="")
hist(data[,2], xlab="Agilent RMA value", breaks=c(0.5*(-20:20)), main="")
myhist(data[,3], xlab="RNAseq RPKM value", maxvalue=100, breakcount=breakcount)
myhist(data[,4], xlab="RNAseq RSEM value", maxvalue=3000, breakcount=breakcount)

glist <- gList()
tmp <- textGrob(label = "a", x = 0.03, y = 0.96, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "b", x = 0.53, y = 0.96, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "c", x = 0.03, y = 0.48, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "d", x = 0.53, y = 0.48, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
grid.draw(glist)

dev.off()
a<-data[,3]
printrange(data[,1], "Affymetrix RMA value");
printrange(data[,2], "Agilent RMA value");
printrange(data[,3], "RNASeq RPKM value");
printrange(data[,4], "RNASeq RSEM value");
