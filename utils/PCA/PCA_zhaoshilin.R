#根据赵世林PCA_2修改。

mypca<-function(filename, number=10,logdata=F,reverse=T,scale=T,col1="red",col2="blue",barplot=F,pos=NULL,cex=1.2,jpeg=F) {
	#读取文件，按Tab分隔，第一列为row name，第一行为col name
	Data <- read.table(filename,header=T,sep="\t",row.names=1);

	If(jpeg == T){
		outfilename = paste(filename,".pca.jpg",sep="");
		jpeg(file=outfilename, width=480, height=640);
	}

	if (barplot==T) {
		op <- par(mfrow=c(2,1), mar=.1+c(2,2,2,2));
	}
	
	if (reverse==T) {
		data<-t(data);
	}
	
	if(logdata == T){
		data<-log(data);
	}
	
	result<-prcomp(data,scale=scale)
	
	plot(result$x[,1],result$x[,2],xlab="PC1",ylab="PC2",main="PCA")
	
	text(result$x[1:number,1],result$x[1:number,2],labels=rownames(result$x)[1:number],col=col1,pos=pos,cex=cex)
	text(result$x[(number+1):length(result$x[,1]),1],result$x[(number+1):length(result$x[,1]),2],labels=rownames(result$x)[(number+1):length(result$x[,1])],col=col2,pos=pos,cex=cex)
	
	abline(h=0,lty=3)
	abline(v=0,lty=3)
	
	if (barplot==T) {
		barplot((result$sdev)^2/sum((result$sdev)^2),names=colnames(result$x),main="Proportion of Variance");
	}
	
	if(jpeg == T){
		dev.off();
	}

	return (result)
}
