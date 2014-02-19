setwd("E:/sqh/Dropbox/papers/published/2013-CNV/paper")

library(grid)

cgh<-read.table("cgh.txt",sep="\t", header=T)
cnvnator<-read.table("cnvnator.txt",sep="\t", header=T)
conifer<-read.table("conifer.txt",sep="\t", header=T)
cnmops<-read.table("cnmops.txt",sep="\t", header=T)

lenarray<-c(nrow(cgh), nrow(cnvnator), nrow(conifer), nrow(cnmops))
maxr = max(lenarray)
out<-matrix(NA, maxr, 4) 
out[1:nrow(cgh),1]<-cgh$length
out[1:nrow(cnvnator),2]<-cnvnator$length
out[1:nrow(conifer),3]<-conifer$length
out[1:nrow(cnmops),4]<-cnmops$length
names<-c("CGH","CNVNator","Conifer","cnmops")
colnames(out)<-names

jpeg(filename="figure1.jpg",width=4000, height=3000,res=300)
par(mfrow = c(1,2))
bplt <- barplot(lenarray, names.arg=names,ylab = "CNV Count", ylim=c(0,maxr*1.1))
text(x=bplt, y=lenarray + 2500, labels=as.character(lenarray))
boxplot(log(out),ylab = "log(CNV length)")

grid.text(label = "(A)", x = 0.01, y = 0.9, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
grid.text(label = "(B)", x = 0.51, y = 0.9, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))

dev.off()

