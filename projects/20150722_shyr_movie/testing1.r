setwd("E:/sqh/programs/r/projects/20150722_shyr_movie/figure")

jpeg("ConceptOfIntervalEstimating_%02d.jpg", width=1000, height=600, res=200)
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

unlink("../ConceptOfIntervalEstimating.mp4")
system("E:/tools/ImageMagick-7.0.0-0-portable-Q16-x64/ffmpeg.exe -r 2 -qscale 2 -i ConceptOfIntervalEstimating_%02d.jpg ../ConceptOfIntervalEstimating.mp4")

jpeg("ConceptOfTestingHypothesis_1_%02d.jpg", width=1000, height=600, res=200)
n<-5;w<-100
x<-matrix( rnorm(n*w, mean=0, sd=1), w, n)
xbar<-round( apply( x, 1, mean),2)
std<-round( apply(x, 1, sd),2 )
z<-round( xbar/(std/sqrt(n)), 2)
alpha.hat<-cumsum( abs(z)>2.78 )
for (i in 1:w){
  stripchart( z[1:i], method = "stack", pch = 16, bty="n", 
              xlab="t", ylim=c(0,10), xlim=c(-6,6), col=4,
              main=paste("data:", paste(round(x[i,],2),collapse=", "),"~ N(0, 1)" ), cex=0.9)
  abline(v=-2.78, col=2, lty=2)
  abline(v=2.78, col=2, lty=2)
  legend("topleft", paste("xbar =", xbar[i]), cex=0.6)
  legend("topright", paste("S =", std[i]),cex=0.6)
  legend("top", paste("t =", z[i] ),cex=0.6 )
  text(-4.5, 7.5, paste("alpha =" , alpha.hat[i], "/", i, "=", round(alpha.hat[i]/i,2) ),cex=0.6 )
}
dev.off()

unlink("../ConceptOfTestingHypothesis_1.mp4")
system("E:/tools/ImageMagick-7.0.0-0-portable-Q16-x64/ffmpeg.exe -r 2 -qscale 2 -i ConceptOfTestingHypothesis_1_%02d.jpg ../ConceptOfTestingHypothesis_1.mp4")

jpeg("ConceptOfTestingHypothesis_2_%02d.jpg", width=1000, height=600, res=200)
x1<-matrix( rnorm(n*w, mean=2, sd=1), w, n)
xbar1<-round( apply( x1, 1, mean),2)
std1<-round( apply(x1, 1, sd),2 )
z1<-round( xbar1/(std1/sqrt(n)), 2)
beta.hat<-cumsum( z1 < 2.78 ) 
for (i in 1:w){
  stripchart( z, method = "stack", pch = 16, bty="n", 
              xlab="t", ylim=c(0,10), xlim=c(-6,15), col=4, 
              main=paste("data:", paste(round(x[i,],2),collapse=", "),"~ N(2, 1)" ), cex=0.9)
  abline(v=-2.78, col=2, lty=2)
  abline(v=2.78, col=2, lty=2)
  stripchart( z1[1:i], method = "stack", pch = 16, 
              bty="n", add=T, col=3, xlim=c(5,15))
  text(10, 8.5, paste("beta =" , beta.hat[i], "/", i, "=", round(beta.hat[i]/i,2) ),cex=0.6 )
  legend("topleft", paste("t =", z1[i] ),cex=0.6 )
}
dev.off()

unlink("../ConceptOfTestingHypothesis_2.mp4")
system("E:/tools/ImageMagick-7.0.0-0-portable-Q16-x64/ffmpeg.exe -r 2 -qscale 2 -i ConceptOfTestingHypothesis_2_%02d.jpg ../ConceptOfTestingHypothesis_2.mp4")
