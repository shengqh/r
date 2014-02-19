# TODO: Add comment
# 
# Author: sqh
###############################################################################


dir="E:/sqh/Science/Project/Shift/10ppm/";

files<-c("Modif.txt","Nomodif.txt","Modif_shift10.txt","Nomodif_shift10.txt");
filelabels<-c("MOD","UNMOD","MOD_SHIFT","UNMOD_SHIFT");

outfilename=paste(dir, "NumMissedCleavage.jpg", sep="");

jpeg(file=outfilename, width=480 * 2, height=640 * 2);

split.screen(c(4, 2));

mybreaks<-c(-100,0.2*-100:100,100);

for(i in 1:4){
  inputfilename = paste(dir, files[i], sep="");
  g<-read.table(inputfilename ,header=T,sep="\t");
  
  for(j in 0:1){
    screennumber <- (i-1) * 2 + (j+1);
    
    screen(screennumber);
    
    title<-paste(filelabels[i]," Charge 2+, NumMissedCleavage ",j,sep="");
    
    t1<-"T-M";
    g1<-subset(subset(subset(subset(g,Decoy=="False"), Modified=="True"), Charge==2), NumMissedCleavage==j);
    t2<-"D-M";
    g2<-subset(subset(subset(subset(g,Decoy=="True"), Modified=="True"), Charge==2), NumMissedCleavage==j);
    t3<-"T-U";
    g3<-subset(subset(subset(subset(g,Decoy=="False"), Modified=="False"), Charge==2), NumMissedCleavage==j);
    t4<-"D-U";
    g4<-subset(subset(subset(subset(g,Decoy=="True"), Modified=="False"), Charge==2), NumMissedCleavage==j);
    
    maxCount<-0;
    if(0 != nrow(g1)){
      r1<-hist(x=g1[,1],breaks=mybreaks,plot=FALSE);
      maxCount<-max(maxCount,r1$counts);
    }
    
    if(0 != nrow(g2)){
      r2<-hist(g2[,1],breaks=mybreaks,plot=FALSE);
      maxCount<-max(maxCount,r2$counts);
    }
    
    if(0 != nrow(g3)){
      r3<-hist(g3[,1],breaks=mybreaks,plot=FALSE);
      maxCount<-max(maxCount,r3$counts);
    }
    
    if(0 != nrow(g4)){
      r4<-hist(g4[,1],breaks=mybreaks,plot=FALSE);
      maxCount<-max(maxCount,r4$counts);
    }
    
    plot(0,xlim=c(-2,10),ylim=c(0,maxCount),type="n",xlab="log(Score)",ylab="Spectrum Count",main=title);
    
    legends<-c();
    colors<-c();
    if(0 != nrow(g1)){
      r1<-hist(g1[,1],breaks=mybreaks,plot=FALSE);
      lines(r1$breaks[2:length(r1$breaks)],r1$counts,col="red");
      legends<-c(legends, t1);
      colors<-c(colors,"red");
    }
    
    if(0 != nrow(g2)){
      r2<-hist(g2[,1],breaks=mybreaks,plot=FALSE);
      lines(r2$breaks[2:length(r2$breaks)],r2$counts,col="blue");
      legends<-c(legends, t2);
      colors<-c(colors,"blue");
    }
    
    if(0 != nrow(g3)){
      r3<-hist(g3[,1],breaks=mybreaks,plot=FALSE);
      lines(r3$breaks[2:length(r3$breaks)],r3$counts,col="green");
      legends<-c(legends, t3);
      colors<-c(colors,"green");
    }
    
    if(0 != nrow(g4)){
      r4<-hist(g4[,1],breaks=mybreaks,plot=FALSE);
      lines(r4$breaks[2:length(r4$breaks)],r4$counts,col="brown");
      legends<-c(legends, t4);
      colors<-c(colors,"brown");
    }
    
    legend(6,maxCount,legend=legends,lty=4,col=colors);
  }
}

close.screen(all = TRUE);
dev.off();
