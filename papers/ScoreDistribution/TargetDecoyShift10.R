# TODO: Add comment
# 
# Author: sqh
###############################################################################


dir="E:/sqh/Science/Project/Shift/10ppm/";

files<-c("Modif.txt","Nomodif.txt","Modif_shift10.txt","Nomodif_shift10.txt");
filelabels<-c("MOD","UNMOD","MOD_SHIFT","UNMOD_SHIFT");
xlabels<-c("log(Score)","-log(ExpectValue)");

outfilename=paste(dir, "shift10.jpg", sep="");

jpeg(file=outfilename, width=480 * 2, height=640 * 2);

split.screen(c(4, 2));

for(i in 1:4){
	inputfilename = paste(dir, files[i], sep="");
	g<-read.table(inputfilename ,header=T,sep="\t");
	
  title<-filelabels[i];
  
  t1<-"T-M";
	g1<-subset(subset(g,Decoy=="False"), Modified=="True");
	t2<-"D-M";
	g2<-subset(subset(g,Decoy=="True"), Modified=="True");
	t3<-"T-U";
	g3<-subset(subset(g,Decoy=="False"), Modified=="False");
	t4<-"D-U";
	g4<-subset(subset(g,Decoy=="True"), Modified=="False");
	
	for(j in 1:2){
		screennumber <- (i-1) * 2 + j;
		
		screen(screennumber);
		
		mybreaks<-c();
		if(1 == j){
			mybreaks<-c(-100,0.2*-100:100,100);
		}
		else{
			mybreaks<-c(-100,0.2*-100:100,100);
		}
		
		maxCount<-0;
		if(0 != nrow(g1)){
			r1<-hist(x=g1[,j],breaks=mybreaks,plot=FALSE);
			maxCount<-max(maxCount,r1$counts);
		}
		
		if(0 != nrow(g2)){
			r2<-hist(g2[,j],breaks=mybreaks,plot=FALSE);
			maxCount<-max(maxCount,r2$counts);
		}
		
		if(0 != nrow(g3)){
			r3<-hist(g3[,j],breaks=mybreaks,plot=FALSE);
			maxCount<-max(maxCount,r3$counts);
		}
		
		if(0 != nrow(g4)){
			r4<-hist(g4[,j],breaks=mybreaks,plot=FALSE);
			maxCount<-max(maxCount,r4$counts);
		}
		
		if(1 == j){
			plot(0,xlim=c(-2,10),ylim=c(0,maxCount),type="n",xlab=xlabels[j],ylab="Spectrum Count",main=title);
		}
		else{
			plot(0,xlim=c(-10,10),ylim=c(0,maxCount),type="n",xlab=xlabels[j],ylab="Spectrum Count",main=title);
		}  
		
		legends<-c();
		colors<-c();
		if(0 != nrow(g1)){
			r1<-hist(g1[,j],breaks=mybreaks,plot=FALSE);
			lines(r1$breaks[2:length(r1$breaks)],r1$counts,col="red");
			legends<-c(legends, t1);
			colors<-c(colors,"red");
		}
		
		if(0 != nrow(g2)){
			r2<-hist(g2[,j],breaks=mybreaks,plot=FALSE);
			lines(r2$breaks[2:length(r2$breaks)],r2$counts,col="blue");
			legends<-c(legends, t2);
			colors<-c(colors,"blue");
		}
		
		if(0 != nrow(g3)){
			r3<-hist(g3[,j],breaks=mybreaks,plot=FALSE);
			lines(r3$breaks[2:length(r3$breaks)],r3$counts,col="green");
			legends<-c(legends, t3);
			colors<-c(colors,"green");
		}
		
		if(0 != nrow(g4)){
			r4<-hist(g4[,j],breaks=mybreaks,plot=FALSE);
			lines(r4$breaks[2:length(r4$breaks)],r4$counts,col="brown");
			legends<-c(legends, t4);
			colors<-c(colors,"brown");
		}
		
		if(1 == j){
			legend(4,maxCount,legend=legends,lty=4,col=colors);
		}
		else{
			legend(0,maxCount,legend=legends,lty=4,col=colors);
		}  
	}
}

close.screen(all = TRUE);
dev.off();
