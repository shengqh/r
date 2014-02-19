estimate_mean<-function(x, col1, col2, binwidth, binsize){
  x <- subset(x, (x[,col1] > 0) & (x[,col2] > 0));
  
  x$logRatio <- log2(x[,col1] / x[,col2]);
  
  f2 <- function(v, c1, c2) {
    ret<-mean(c(v[c1], v[c2]));
    return(ret);
  }
  x$aveIntensity <- apply(x, 1, f2, c1=col1, c2=col2);
  x$binbox <- 1;
  
  minintensity = min(x$aveIntensity);
  maxintensity = max(x$aveIntensity);
  binindex<-data.frame();
  binhigh<-maxintensity;
  while(binhigh > minintensity){
    rs <- data.frame();
    binlow<-binhigh;
    while(nrow(rs) < binsize){
      binlow<-binlow-binwidth;
      rsleft<-subset(x, aveIntensity <= binlow);
      
      if(nrow(rsleft) < binsize){
        binlow<-minintensity - 0.01;
      }
      
      rs <- subset(x, (aveIntensity <= binhigh) & (aveIntensity > binlow));
      
      if(binlow < minintensity){
        break;
      }
    }
    
    rsv <- as.vector(as.matrix(rs$logRatio));
    rsv99 <- get_data_99(rsv);
    binmean <- mean(rsv99);
    binstdev <-get_stdev_with_mean(rsv99, binmean);
    #defaultstdev<-sd(rsv99)
    
    #if(binstdev != defaultstdev){
    #cat("binstdev=", binstdev, " ; defaultstdev=", defaultstdev, "\n")
    #}
    
    binindex<-rbind(c(binlow, binhigh, binmean, binstdev), binindex);
    
    x$binbox[(x$aveIntensity <= binhigh) & (x$aveIntensity > binlow)] <- nrow(binindex)
    binhigh<-binlow;
  }
  
  x$binbox <- nrow(binindex) + 1 - x$binbox
  
  colnames(binindex) <- c("binlow", "binhigh", "binmean", "binstdev")
  
  return (list(value = x, bins=binindex))
}

#get middle 99% data
get_data_99<-function(rsv){
  sortedrsv <- sort(rsv)
  
  ss<-max(2, trunc(length(rsv) * 0.005))
  se<-trunc(length(rsv) * 0.995)
  
  rsv99<-sortedrsv[ss:se]
  
  return (rsv99)
}

#get stdev with assigned mean
get_stdev_with_mean<-function(rdata, rmean){
  stdev <- sqrt( sum((rdata - rmean) ** 2 ) / (length(rdata) - 1) )
  return (stdev)
}

calc_stdev<-function(bin, xx){
  binlow<-as.numeric(bin["binlow"])
  binhigh<-as.numeric(bin["binhigh"])
  allmean<-as.numeric(bin["mean"])
  
  rs <- subset(xx, (aveIntensity <= binhigh) & (aveIntensity > binlow))
  
  rsv <- as.vector(as.matrix(rs$logRatio))
  
  rsv99 <- get_data_99(rsv)
  
  stdev<-get_stdev_with_mean(rsv99, allmean)
  
  return(stdev)
}

calc_probability<-function(v, bins){
  logratio<-v["logRatio"]
  binbox<-v["binbox"]
  
  bmean<-bins[binbox,"mean"]
  bstdev<-bins[binbox,"stdev"]
  
  plower<-pnorm(logratio, mean=bmean, sd=bstdev, lower.tail=T);
  phigher<-pnorm(logratio, mean=bmean, sd=bstdev, lower.tail=F);
  probability<- (min(plower, phigher) * 2);
  
  return (probability)
}

setwd("D:\\sqh\\programs\\r\\Guoyan\\Shanks")
filename<-"HNSC_RSeqV2.20121119.tsv";
x=read.table(filename, sep="\t", header=T, row.names=1)
x$tail=apply(x[,2:10], 1, max)

col1<-1;
col2<-12;
  
binwidth=0.05
binsize=50
  
ret<-estimate_mean(x, col1, col2, binwidth, binsize);
  
xx<-ret$value
bins<-ret$bins
  
bins$mean <-mean(bins$binmean)
  
bins$stdev <- apply(bins, 1, calc_stdev, xx)
  
xx$probability<-apply(xx, 1, calc_probability, bins)
  
  png(filename=paste0(filename, ".", binwidth, ".", binsize, ".png"),width=4000,height=3000,units="px",res=300);
  cols<-rep("black",nrow(xx));
  cols[xx$probability <= 0.05]<-'red';
  cols[xx$probability <= 0.01]<-'blue';
  plot(xx$aveIntensity, xx$logRatio, pch='+', col=cols, xlab="Average Normalized Density", ylab="log(SH1/SH2)");
  dev.off();
  
  write.csv(xx, file=paste0(filename, ".", binwidth, ".", binsize, ".csv"));
  write.csv(bins, file=paste0(filename, ".", binwidth, ".", binsize, ".bin.csv"));
