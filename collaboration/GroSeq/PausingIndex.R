setwd("G:/sqh/projects/GRO-Seq/data/figure")

readPausingData<-function(filename, breaks, pvalue){
  data<-read.csv(filename);
  
  all<-subset(data, data[,2] < 100)
  all<-subset(all, data[,2] > 0)
  alldata<-all[,2]
  allhist<-hist(alldata,breaks=breaks,plot=FALSE)
  
  sig<-subset(data, data[,3] < pvalue)
  sig<-subset(sig, data[,2] < 100)
  sig<-subset(sig, data[,2] > 0)
  sigdata<-sig[,2]
  sighist<-hist(sigdata,breaks=breaks,plot=FALSE)
  
  result<-list(alldata=alldata, allhist=allhist, sigdata=sigdata, sighist=sighist)
  return (result)
}

breaks=c(0:100)

wtdata<-readPausingData("WT.csv", breaks, 0.01)
nulldata<-readPausingData("Null.csv", breaks, 0.01)

maxCount<-max(wtdata$allhist$counts, wtdata$sighist$counts,nulldata$allhist$counts,nulldata$sighist$counts);

plot(0,xlim=c(0,100),ylim=c(0,maxCount),type="n",xlab="Pausing Index",ylab="Number of genes",main="Pausing Index Distribution");

hist(wtdata$alldata, breaks=breaks, col="grey", add=TRUE);
hist(wtdata$sigdata, breaks=breaks, col="red", add=TRUE, legend="WT-Significantly paused (p val < 0.01)");
lines(nulldata$allhist$breaks[2:length(nulldata$allhist$breaks)],nulldata$allhist$counts,col="blue");
lines(nulldata$sighist$breaks[2:length(nulldata$sighist$breaks)],nulldata$sighist$counts,col="green");




