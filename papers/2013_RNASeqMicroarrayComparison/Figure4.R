library(grid)

rootdir<-"E:/sqh/Dropbox/papers/2013-RNAseqMicroarrayComparison";
datafile<-paste0(rootdir, "/Figure3_BRCA_PValue.csv");
data<-read.csv(datafile, header=T, row.names=1, check.names=F)
data<-subset(data, !is.na(data[,6]))
sunion<-data[order(data$RPKMLog2Mean),]
step=nrow(sunion) / 4
step1<-sunion[1:(step),]
step2<-sunion[(step+1):(2*step),]
step3<-sunion[(2*step+1):(3*step),]
step4<-sunion[(3*step+1):nrow(sunion),]

tiff(file=paste0(rootdir,"/Figure4_IntensityRelatedFoldChange.tif"),width=2000, height=2000, res=300, compression="lzw")
par(mar=c(6,5,2,1))
par(mfrow=c(2,2))

drawplot<-function(d, xlab){
  plot(d[,"RPKMFold"], d[,"AgilentFold"], xlab=xlab, ylab="log2(Agilent fold change)", xlim=c(-10,10), ylim=c(-10,10))
  fit = lm(AgilentFold ~ RPKMFold, data=d)
  fitsum = summary(fit)
  r2 = fitsum$adj.r.squared
  lm_coef <- coef(fit)
  rp = vector('expression',2)
  if (lm_coef[1] > 0){
    rp[1] = substitute(expression(y == VALUE1*x + VALUE2), 
                     list(VALUE1 = format(lm_coef[2],dig=3),
                          VALUE2 = format(lm_coef[1],dig=3)
                          ))[2]
  }
  else{
    rp[1] = substitute(expression(y == VALUE1*x - VALUE2), 
                       list(VALUE1 = format(lm_coef[2],dig=3),
                            VALUE2 = format(abs(lm_coef[1]),dig=3)
                       ))[2]
  }
  rp[2] = substitute(expression(italic(R)^2 == MYVALUE), 
                     list(MYVALUE = format(r2,dig=3)))[2]
  
  abline(fit)
  legend("topleft", bty="n", legend=rp)
}

drawplot(step1, "log2(RPKM[0-25%] fold change)")
drawplot(step2, "log2(RPKM[25-50%] fold change)")
drawplot(step3, "log2(RPKM[50-75%] fold change)")
drawplot(step4, "log2(RPKM[75-100%] fold change)")

glist <- gList()
tmp <- textGrob(label = "a", x = 0.03, y = 0.96, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "b", x = 0.53, y = 0.96, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "c", x = 0.03, y = 0.46, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "d", x = 0.53, y = 0.46, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
grid.draw(glist)


dev.off();