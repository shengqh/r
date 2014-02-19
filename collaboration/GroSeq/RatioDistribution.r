pvaluecalc<-function(v,binsize,x){
  ranki<-v[ncol(x)];
  ratioi<-v[ncol(x)-2];
  binlow<-ranki-binsize;
  binhigh<-ranki+binsize;
  rs <- subset(x, (x["rankAveIntensity"] >= binlow));
  rs <- subset(rs, (rs["rankAveIntensity"] <= binhigh));
  rsv <- as.vector(as.matrix(rs["logRatio"]));
  count <- length(rsv);
  rave <- mean(rsv);
  rstdev <- sd(rsv);
  if(rstdev == 0){
    pvalue=1;
  }
  else{
    pvaluelower<-pnorm(ratioi, mean=rave, sd=rstdev, lower.tail=T);
    pvaluehigher<-pnorm(ratioi, mean=rave, sd=rstdev, lower.tail=F);
    pvalue<- (min(pvaluelower, pvaluehigher) * 2);
  }
  return(c(count, rave, rstdev, pvalue));
}

pvaluecalc_byabundance<-function(v,binwidth,binsize,x){
  ratioi<-v['logRatio'];

  aveIntensity<-v['aveIntensity'];
  minIntensity<-aveIntensity - binwidth;
  maxIntensity<-aveIntensity + binwidth;
  rs <- subset(x, (x["aveIntensity"] >= minIntensity));
  rs <- subset(rs, (rs["aveIntensity"] <= maxIntensity));

  if(nrow(rs) < binsize){
    ranki<-v['rankAveIntensity'];

    binlow<-ranki-binsize/2;
    binhigh<-ranki+binsize/2;
    rs <- subset(x, (x["rankAveIntensity"] >= binlow));
    rs <- subset(rs, (rs["rankAveIntensity"] <= binhigh));
  }
  
  rsv <- as.vector(as.matrix(rs["logRatio"]));
  count <- length(rsv);
  rave <- mean(rsv);
  rstdev <- sd(rsv);
  if(rstdev == 0){
    pvalue=1;
  }
  else{
    pvaluelower<-pnorm(ratioi, mean=rave, sd=rstdev, lower.tail=T);
    pvaluehigher<-pnorm(ratioi, mean=rave, sd=rstdev, lower.tail=F);
    pvalue<- (min(pvaluelower, pvaluehigher) * 2);
  }
  return(c(count, rave, rstdev, pvalue));
}

abundance_different <-function(x, col1, col2, binwidth, binsize){
  x1<- (x[,col1] > 0);
  x <- subset(x, x1);

  x2<- (x[,col2] > 0);
  x <- subset(x, x2);
  
  z <- log(x[,col1] / x[,col2]);
  x <- cbind(x, z);
  colnames(x)[ncol(x)] <- "logRatio";
  
  f2 <- function(v, c1, c2) {
    ret<-log(mean(c(v[c1], v[c2])));
    return(ret);
  }
  a <- apply(x, 1, f2, c1=col1, c2=col2);
  x <- cbind(x, a);
  colnames(x)[ncol(x)] <- "aveIntensity";
  
  r<-rank(a)
  x <- cbind(x, r);
  colnames(x)[ncol(x)] <- "rankAveIntensity";
  
  p<-t(apply(x,1,pvaluecalc_byabundance,binwidth,binsize,x));
  x<-cbind(x,p);
  colnames(x)[ncol(x)-3] <- "count";
  colnames(x)[ncol(x)-2] <- "mean";
  colnames(x)[ncol(x)-1] <- "stdev";
  colnames(x)[ncol(x)] <- "pvalue";

  q<-p.adjust(p[,4])
  x<-cbind(x,q);
  colnames(x)[ncol(x)] <- "fdr";
  return (x);
}

for(filename in c(
  "H:/shengquanhu/projects/GRO-Seq/20121101/groseq-HDAC3-uniq.rRNA.gb500bp-promoter50-pp-called-forstat",
  "H:/shengquanhu/projects/GRO-Seq/20121101/groseq-HDAC3-uniq.rRNA.gb500bp-promoter50-pp-allrefseqs-forstat")){
  x<-read.csv(paste0(filename,".csv"),header=T,row.names=1)

  col1<-4;
  col2<-9;

  binwidth=0.05
  binsize=100
  #for(binwidth in c(0.05,0.1,0.2)){
    #for (binsize in c(30,50,100,150,200)){
      xx<-abundance_different(x, col1, col2, binwidth, binsize);
      png(filename=paste0(filename, ".", binwidth, ".", binsize, ".png"),width=4000,height=3000,units="px",res=300);
      cols<-rep("black",nrow(xx));
      cols[xx['pvalue'] <= 0.05]<-'red';
      cols[xx['pvalue'] <= 0.01]<-'blue';
      plot(xx[,12], xx[,11], pch='+', col=cols, xlab="Average Normalized Density", ylab="log(SH1/SH2)");
      dev.off();
  
      write.csv(xx, file=paste0(filename, ".", binwidth, ".", binsize, ".csv"));
    #}
  #}
}