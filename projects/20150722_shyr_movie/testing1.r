setwd("E:/sqh/programs/r/projects/20150722_shyr_movie/figure")
jpeg("figure1_%02d.jpg", width=1000, height=600, res=200)
par(mfcol=c(1,2), mar = c(3, 2, 0.5, 1))
n<-30;w<-100
alpha<-0.05
set.seed(1235)
x<-matrix( rnorm(n*w), w, n)
ci<-function(a) cbind(mean(a)+qnorm(alpha/2)*sqrt(1/n), mean(a)-qnorm(alpha/2)*sqrt(1/n)) 
out<-t(apply(x, 1, ci))
for (i in 1:w) {
  plot(0.5-dnorm(seq(-3,3,0.1)), seq(-3,3,0.1), type="l", bty="n", 
       yaxt="n", ylab="", xlab="", xaxt="n" )
  axis(4)
  points(0.506, 0, pch=15, col=4)
  plot(rep(1,2), out[1,], xlim=c(1, 100), type="n", 
       ylim=c(-1.0, 1.0), xlab="", ylab="", bty="n", yaxt="n")
  out1<-cbind( out, seq(1, w), seq(1,w))
  for (j in 1:i){
    ifcol<-ifelse( out[j,1]<0 & out[j,2]>0, 1, 2)
    arrows(out1[j,3], out1[j,1], out1[j,4], out1[j,2], code=3, angle=90, length=0.02, col=ifcol)
  }
  abline(h=0, lty=2, col="blue", lwd=2)
}
dev.off()

unlink("figure1.mp4")
system("E:/tools/ImageMagick-7.0.0-0-portable-Q16-x64/ffmpeg.exe -r 2 -qscale 2 -i figure1_%02d.jpg figure1_2pics_per_second.mp4")
