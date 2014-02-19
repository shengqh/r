# TODO: Add comment
# 
# Author: sqh
###############################################################################

#dir="E:/sqh/Science/Project/Shift/sequest_10ppm/";

#inputfilename=paste(dir,"Modif.txt",sep="");

dir="E:/sqh/Science/Project/wyb/";
inputfilename=paste(dir,"3t3l1_phospho.txt",sep="");

outfilename=paste(dir,"Score.jpg",sep="");

g<-read.table(inputfilename ,header=T,sep="\t");

jpeg(file=outfilename, width=480 * 3, height=640 * 3);

split.screen(c(3, 3));

mybreaks<-c(-100,0.2*-40:40,100);

screenindex = 0;
charge = 2;
miss = 0;

    screenindex<-screenindex+1;
    screen(screenindex);
    
    gChargeMiss<-subset(subset(g, Charge==charge), NumMissedCleavage==miss);
    title<-paste("Charge",charge,"; NumMissCleavage",miss,"; Count",nrow(gChargeMiss));
    
    gt<-subset(gChargeMiss,Decoy=="False");
    
    gt1<-subset(gt, Modified=="True");
    rt1<-hist(gt1[,1],breaks=mybreaks,plot=FALSE);
    
    gt2<-subset(gt, Modified=="False");
    rt2<-hist(gt2[,1],breaks=mybreaks,plot=FALSE);
    
    gd<-subset(gChargeMiss,Decoy=="True");
    
    gd1<-subset(gd, Modified=="True");
    rd1<-hist(gd1[,1],breaks=mybreaks,plot=FALSE);
    
    gd2<-subset(gd, Modified=="False");
    rd2<-hist(gd2[,1],breaks=mybreaks,plot=FALSE);
    
    maxCount<-max(rt1$counts, rt2$counts, rd1$counts, rd2$counts);
    
    plot(0,xlim=c(-8,8),ylim=c(0,maxCount),type="n",xlab="log(Score)",ylab="Spectrum Count",main=title);
    
    lines(rt1$breaks[2:length(rt1$breaks)],rt1$counts,col="red");
    lines(rt2$breaks[2:length(rt2$breaks)],rt2$counts,col="blue");
    lines(rd1$breaks[2:length(rd1$breaks)],rd1$counts,col="green");
    lines(rd2$breaks[2:length(rd2$breaks)],rd2$counts,col="brown");
    
    legend(-8,maxCount,legend=c(
            paste("MOD_T",nrow(gt1)),
            paste("MOD_D",nrow(gd1)),
            paste("UNMOD_T",nrow(gt2)),
            paste("UNMOD_D",nrow(gd2))
            ),lty=4,col=c("red","green","blue","brown"));

close.screen(all = TRUE);
dev.off();
