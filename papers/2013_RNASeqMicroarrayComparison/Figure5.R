library(grid)

drawplot<-function(d, xlab){
  plot(d[,"RPKM"], d[,"RSEM"], xlab=xlab, ylab="log2(RSEM)"
       , xlim=c(-18,18), ylim=c(-18,18)
  )
  par(xpd=FALSE)
  fit = lm(RPKM ~ RSEM , data=d)
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

rootdir<-"H:/sqh/Dropbox/Projects/TCGA_Comparison";
datafile<-paste0(rootdir, "/Figure5_Exon_RPKM_RSEM.median.tsv");
data<-read.table(datafile, header=T, row.names=1, check.names=F, sep="\t")

cat("Median of exon length = ", median(data$Length))
cat("Mean of exon length = ", mean(data$Length))

data$Length<-log10(data$Length)

tiff(file=paste0(rootdir,"/Figure5_ExonCorrelation.tif"),width=3000, height=2000, res=300, compression="lzw")
par(mar=c(6,5,2,1))
par(mfrow=c(2,3))

##hist of length
hh<-hist(data$Length, xlab="log10(Exon length)", main="", breaks=20)

breaks=hh$breaks
rsem0<-data[(data$RSEM==0) & (data$RPKM != 0),]
rpkm0<-data[(data$RPKM==0) & (data$RSEM !=0),]
rsem0count<-c()
rpkm0count<-c()

for(i in 1:(length(breaks)-1)){
  c1<-rsem0[(rsem0$Length >= breaks[i]) & (rsem0$Length < breaks[i+1]),]
  rsem0count<-c(rsem0count, nrow(c1)) 
  c2<-rpkm0[(rpkm0$Length >= breaks[i]) & (rpkm0$Length < breaks[i+1]),]
  rpkm0count<-c(rpkm0count, nrow(c2)) 
}

##rsem/rpkm differential
ymax<-max(rsem0count, rpkm0count) * 1.2
plot(0,bty='n',pch='',ylab='Frequency',xlab='log10(Exon length)',xlim=c(min(breaks), max(breaks)),ylim=c(0, ymax))
lines(x=breaks[1:(length(breaks)-1)], y=rsem0count, col="red")
lines(x=breaks[1:(length(breaks)-1)], y=rpkm0count, col="blue")
legend("topleft", legend=c("RPKM(RSEM==0)", "RSEM(RPKM==0)"), text.col=c("red","blue"), bty='n')

nozero<-data[(data$RSEM != 0) & (data$RPKM != 0),]
nozero$RPKM<-log2(nozero$RPKM)
nozero$RSEM<-log2(nozero$RSEM)
#drawplot(nozero, "log2(RPKM)", "c")

lbreaks=c(0, breaks[6:18], breaks[length(breaks)])
r2list<-c()
x<-c()
for(i in 1:(length(lbreaks)-1)){
  dd<-nozero[(nozero$Length >= lbreaks[i]) & (nozero$Length < lbreaks[i+1]),]
  fit = lm(RPKM ~ RSEM, data=dd)
  fitsum = summary(fit)
  r2 = fitsum$adj.r.squared
  r2list<-c(r2list,r2)
  x<-c(x, mean(c(lbreaks[i],lbreaks[i+1])))
}
##R2 related to exon length
plot(0,bty='n',pch='',ylab='R2 (RSEM ~ RPKM)',xlab='log10(Exon length)',xlim=c(min(lbreaks), max(lbreaks)),ylim=c(0, 1))
lines(x=x, y=r2list, col="black", type='o')


low<-nozero[(nozero$Length >= 0) & (nozero$Length <= log10(20)),]
drawplot(low,paste0("log2(RPKM) [exon length <= 20, N=", nrow(low), "]"))
low<-nozero[(nozero$Length > log10(20)) & (nozero$Length <= log10(50)),]
drawplot(low,paste0("log2(RPKM) [exon length=20~50, N=", nrow(low), "]"))
low<-nozero[(nozero$Length > log10(50)),]
drawplot(low,paste0("log2(RPKM) [exon length >= 50, N=", nrow(low), "]"))

glist <- gList()
tmp <- textGrob(label = "a", x = 0.015, y = 0.96, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "b", x = 0.345, y = 0.96, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "c", x = 0.68, y = 0.96, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "d", x = 0.015, y = 0.48, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "e", x = 0.345, y = 0.48, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)
tmp <- textGrob(label = "f", x = 0.68, y = 0.48, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
glist <- gList(glist, tmp)

grid.draw(glist)


dev.off();