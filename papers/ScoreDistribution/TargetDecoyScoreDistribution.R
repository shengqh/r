# TODO: Add comment
# 
# Author: sqh
###############################################################################


#inputfilename="E:/sqh/Science/Project/Shift/10ppm/modif.txt";
#inputfilename="E:/sqh/Science/Project/Shift/sequest_10ppm/modif.txt";
#inputfilenames=c("Z:/Mouse_Liver/new_build/summary_0.5_0.1/O_SCX_buildsummary_FDR0.5_0.1/Cell/SCX_O_Cell.unique1.unique.txt",
#    "Z:/Mouse_Liver/new_build/summary_0.5_0.1/O_SCX_buildsummary_FDR0.5_0.1/Cell/SCX_O_Cell.unique1.all.txt");
inputfilenames=c("E:/sqh/Science/Project/wyb/3t3l1_phospho.txt");

for(inputfilename in inputfilenames){
outfilename=paste(inputfilename,".jpg",sep="");

g<-read.table(inputfilename ,header=T,sep="\t");

jpeg(file=outfilename, width=480 * 3, height=640 * 3);

split.screen(c(3, 3));

mybreaks<-c(-100,0.2*-40:40,100);

screenindex = 0;
for(charge in 1:3){
  for(miss in 0:2){
    screenindex<-screenindex+1;
    screen(screenindex);
    
    gChargeMiss<-subset(subset(g, Charge==charge), NumMissedCleavage==miss);
    if(nrow(gChargeMiss) > 0){
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
    
    legends<-c();
    colors<-c();
    
    if(nrow(gt1) > 0){
      lines(rt1$breaks[2:length(rt1$breaks)],rt1$counts,col="red");
      legends<-c(legends,paste("MOD_T",nrow(gt1)));
      colors<-c(colors,"red");
    }
    
    if(nrow(gd1) > 0){
      lines(rd1$breaks[2:length(rd1$breaks)],rd1$counts,col="blue");
      legends<-c(legends,paste("MOD_D",nrow(gd1)));
      colors<-c(colors,"blue");
    }
    
    if(nrow(gt2) > 0){
      lines(rt2$breaks[2:length(rt2$breaks)],rt2$counts,col="green");
      legends<-c(legends,paste("UNMOD_T",nrow(gt2)));
      colors<-c(colors,"green");
    }
    
    if(nrow(gd2) > 0){
      lines(rd2$breaks[2:length(rd2$breaks)],rd2$counts,col="brown");
      legends<-c(legends,paste("UNMOD_D",nrow(gd2)));
      colors<-c(colors,"brown");
    }
    
      legend(-8,maxCount,legends ,lty=4,col=colors);
    }
  }
}

close.screen(all = TRUE);
dev.off();
}