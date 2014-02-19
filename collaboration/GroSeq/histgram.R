myhist<-function(source_data, col_names, col_index, target_dir, breaks, logneeded, newpng){
	mydata<-na.omit(source_data[,col_index])
  
  mydata<-mydata[mydata > 0]
  
	if (logneeded) {
	  mydata<-log(mydata)
    lab=paste0("log(", col_names[col_index], ")")
	}
  else{
    lab=col_names[col_index]
  }
  
	
	if (newpng){
    if (logneeded){
      targetfile<-paste0(target_dir, col_names[col_index], ".log.png")
    }
    else{
      targetfile<-paste0(target_dir, col_names[col_index], ".png")
    }
    png(filename=targetfile, width=1024, height=768)
	}
  
  mydata <- mydata[mydata > min(breaks)]
	mydata <- mydata[mydata < max(breaks)]
	
  h<-hist(mydata, breaks, xlab=lab, main=paste0("Histogram of ",col_names[col_index]), freq=TRUE)
  
  if (logneeded){
    xfit<-seq(min(mydata), max(mydata),length=100)
    yfit<-dnorm(xfit, mean=mean(mydata), sd=sd(mydata))
    yfit<-yfit*diff(h$mids[1:2])*length(mydata)
    lines(xfit,yfit,col="blue",lwd=2)
  }

	if (newpng){
	      dev.off()
	}
}

mydir<-"D:/sqh/projects/GRO-Seq/data/"

#sourcefile<-paste0(mydir,"MTG8-density-forstatistics.csv")
sourcefile<-paste0(mydir,"MTG8-gro-uniq-mappable-count-forstat.csv")

c<-read.table(sourcefile,header=T,sep=",")

mylogbreaks<-c(0.2 * -25:75)

mybreaks<-c(-10:300)

png(filename=paste0(sourcefile,".png"), width=800 * 2, height=600 * 2)

split.screen(c(2,2))

screen(1)

myhist(c, colnames(c), 4, mydir, mylogbreaks, 1, 0)

screen(2)

myhist(c, colnames(c), 7, mydir, mylogbreaks, 1, 0)

screen(3)

myhist(c, colnames(c), 4, mydir, mybreaks, 0, 0)

screen(4)

myhist(c, colnames(c), 7, mydir, mybreaks, 0, 0)

close.screen(all = TRUE)

dev.off()
