  getPausingData<-function(data, breaks, pvalue){
    wtall<-subset(data, data[,10] < 100)
    wtall<-subset(wtall, wtall[,10] > 0)
    wtalldata<-wtall[,10]
    wtallhist<-hist(wtalldata,breaks=breaks,plot=FALSE)
    
    wtsig<-subset(data, data[,20] < pvalue)
    wtsig<-subset(wtsig, wtsig[,10] < 100)
    wtsig<-subset(wtsig, wtsig[,10] > 0)
    wtsigdata<-wtsig[,10]
    wtsighist<-hist(wtsigdata,breaks=breaks,plot=FALSE)
    
    wtdata<-list(alldata=wtalldata,allhist=wtallhist,sigdata=wtsigdata,sighist=wtsighist);
    
    nullall<-subset(data, data[,19] < 100)
    nullall<-subset(nullall, nullall[,19] > 0)
    nullalldata<-nullall[,19]
    nullallhist<-hist(nullalldata,breaks=breaks,plot=FALSE)
    
    nullsig<-subset(data, data[,21] < pvalue)
    nullsig<-subset(nullsig, nullsig[,19] < 100)
    nullsig<-subset(nullsig, nullsig[,19] > 0)
    nullsigdata<-nullsig[,19]
    nullsighist<-hist(nullsigdata,breaks=breaks,plot=FALSE)
    
    nulldata<-list(alldata=nullalldata,allhist=nullallhist,sigdata=nullsigdata,sighist=nullsighist);
    
    result<-list(wtdata=wtdata, nulldata=nulldata)
    
    return (result)
  }

  mydir<-"D:/Sync/Dropbox/GroSeq/20120907/"

  setwd(mydir)

files<-c("MTG8-gro-uniq-mapped_promoter50bp-count-fdr.csv","HDAC3-gro-uniq-mapped_promoter50bp-count-fdr.csv")

names<-c("WT","CTRL")

dsnames<-c("MTF8","HDAC3")

for(f in 1:length(files)){
  data<-read.csv(files[f])
  
  breaks=c(0:100)
  
  fdata<-getPausingData(data, breaks, 0.01)
  
  wtdata<-fdata$wtdata;
  nulldata<-fdata$nulldata;
  
  maxCount<-max(wtdata$allhist$counts, 
                nulldata$allhist$counts)
  
  overlap = wtdata$sighist;
  for(i in 1:length(overlap$counts)){ 
    if(wtdata$sighist$counts[i] > 0 & nulldata$sighist$counts[i] > 0){
      overlap$counts[i] = min(wtdata$sighist$counts[i],nulldata$sighist$counts[i])
    } else {
      overlap$counts[i] = 0
    }
  }
  
  colors<-c("black","red","blue", "green", rgb(0.9,0.9,1.0))
  
  tiff(filename=paste0(files[f],".tiff"),width=3000,height=2000,compression="lzw",res=300)
  plot(0,xlim=c(0,100),ylim=c(0,maxCount),type="n",xlab="Pausing Index",ylab="Number of genes",main=paste0(dsnames[f], " Pausing Index Distribution"));
  plot(wtdata$sighist,col=colors[2],add=TRUE)
  plot(nulldata$sighist,col=colors[4],add=TRUE)
  plot(overlap,col=colors[5],add=TRUE)
  lines(wtdata$allhist$breaks[2:length(wtdata$allhist$breaks)],wtdata$allhist$counts,col=colors[1],lwd=2)
  lines(nulldata$allhist$breaks[2:length(nulldata$allhist$breaks)],nulldata$allhist$counts,col=colors[3],lwd=2)
  
  legends<-c(paste0(names[f], "-Allgenes"),paste0(names[f],"-Significantly paused (fdr < 0.05)"),
		  "Null-Allgenes","Null-Significantly paused (fdr < 0.05)",
		  paste0("Overlap count between ", names[f], " and Null"))
  
  legend("topright",maxCount,legend=legends,lty=1,col=colors);
  dev.off()
}
