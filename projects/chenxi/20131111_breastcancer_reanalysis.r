########################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("genefilter")
#biocLite("hgu133a.db")
#install.packages("sfsmisc")
########################################################

library("affy");
library("preprocessCore");
library("sfsmisc")
library("genefilter")
library("hgu133a.db")

#########################################################
## function corresponds to function gap in package SAGx
#########################################################
gapStat <- function (data, class = rep(1, nrow(data)), M = 500){
  if (!(length(class) == nrow(data)))
    stop("Length of class vector differs from nrow of data")
  if(M <= 0)
    stop("'M' has to be a positive integer")
  
  data <- as.matrix(data)
  data <- scale(data, center = TRUE, scale = FALSE)
  M <- trunc(M)
  
  pw.dist <- function(x){ sum(dist(x)/ncol(x))/2 }
  
  temp1 <- log(sum(by(data, factor(class), pw.dist)))
  veigen <- svd(data)$v
  x1 <- crossprod(t(data), veigen)
  z1 <- matrix(data = NA, nrow = nrow(x1), ncol = ncol(x1))
  
  tots <- vector(length = M)
  for (k in 1:M) {
    for (j in 1:ncol(x1)) {
      min.x <- min(x1[, j])
      max.x <- max(x1[, j])
      z1[, j] <- runif(nrow(x1), min = min.x, max = max.x)
    }
    z <- crossprod(t(z1), t(veigen))
    tots[k] <- log(sum(by(z, factor(class), pw.dist)))
  }
  out <- c(mean(tots) - temp1, sqrt(1 + 1/M) * sd(tots))
  names(out) <- c("Gap statistic", "SE of simulation")
  
  return(out)
}

#########################################################
## clustering via function kmeans of package stats where number of centers is 
## choosen via the gap statistic
#########################################################
## arguments identical to function kmeans but without argument "centers"
## B: integer, number of Monte Carlo samples
## hclust: logical choose centers via function hclust using method "average" (cf. 
## p. 318 in VR (2002))
kmeansGap <- function(x, iter.max = 10, nstart = 1, 
                      algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                      k.max = 20, M = 100, hclust = FALSE){
  k.max <- trunc(k.max)
  if(k.max < 2) stop("'k.max' has to be >= 2")
  
  k <- 1
  km.new <- NULL
  gap.new <- gapStat(data = x, class = rep(1, nrow(x)), M = M)
  gap.stat <- gap.new
  
  repeat{
    km.old <- km.new
    gap.old <- gap.new
    k <- k + 1
    
    if(k > k.max){
      warning("'k.max' reached kmeans result for k.max centers returned.")
      break
    }
    
    if(hclust){
      hc <- hclust(dist(x), method = "average")
      initial <- tapply(x, list(rep(cutree(hc, k), ncol(x)), col(x)), mean)
      dimnames(initial) <- list(NULL, dimnames(x)[[2]])
      km.new <- kmeans(x, centers = initial, iter.max = iter.max, nstart = nstart, 
                       algorithm = algorithm)
    }else{
      km.new <- kmeans(x, centers = k, iter.max = iter.max, nstart = nstart, 
                       algorithm = algorithm)
    }
    gap.new <- gapStat(data = x, class = km.new$cluster, M = M)
    gap.stat <- rbind(gap.stat, gap.new)
    
    if(gap.old[1] - gap.new[1] + gap.new[2] >= 0)
      break
  }
  rownames(gap.stat) <- NULL
  gap.stat <- cbind(1:nrow(gap.stat), gap.stat)
  colnames(gap.stat) <- c("number of clusters", "gap statistic", "SE of simulation")
  
  if(k == 2){ 
    warning("number of clusters equal to 1 => NULL is returned")
    return(NULL)
  }
  
  res <- list(gapStat = gap.stat, cluster = km.old)
  class(res) <- "clusterGap"
  
  return(res)
}

plot.clusterGap <- function(x, ...){
  plot(1:nrow(x$gapStat), x$gapStat[,2], ylab = "gap statistic", xlab = "number of clusters", type = "l")
  plotrix::plotCI(1:nrow(x$gapStat), x$gapStat[,2], x$gapStat[,3], add = TRUE)
  title("Gap statistic for clustering")
}

########main procedure########

if (.Platform$OS.type == "unix") {
  setwd("/scratch/cqs/shengq1/breastcancer/final/trainingdataset")
}else{
  setwd("H:/shengquanhu/projects/BreastCancer/final/20131111_brca_microarray_reanalysis");
}

datafiles<-c("breastcancer_affymetrix_ar711.tsv", 
             "breastcancer_affymetrix_1058.tsv"
);
for(datafile in datafiles){
  qnormfile<-paste0(datafile,".qnorm");
  bnormfile<-paste0(qnormfile,".bnorm.genes")
  sdfile<-paste0(bnormfile,".sd0.8")
  sd_rdata<-paste0(sdfile, ".RData")
  kmeans_rdata<-paste0(sdfile, ".kmeans.RData")
  if(file.exists(kmeans_rdata)){
    load(kmeans_rdata)
  }else{
    if(!file.exists(sd_rdata)){
      cat(paste0("reading file ", datafile, " ...\n"))
      
      #read/transfer data
      data<-read.table(datafile,header=T,row.name=1);
      tdata<-t(data);
      rm(data);
      
      cat("quantile normalization ...\n")
      
      normalize.quantiles.robust(x=tdata,copy=F);
      qnormdata<-log2(tdata);
      rm(tdata);
      
      save(qnormdata,file=paste0(qnormfile,".RData"));
      
      cat("batch normalization by linear model ...\n")
      
      batchdefinitionfile<-paste0(datafile,".design");
      batch<-read.table(batchdefinitionfile)
      foo<-lm(t(qnormdata) ~ as.factor(batch$V1))
      rm(qnormdata)
      norm_all<-t(foo$res)
      rm(foo)
      nm <- as.vector(sapply(rownames(norm_all), function(x) substr(x,3,nchar(x)),simplify=TRUE))
      rownames(norm_all)<-nm
      
      cat("filter gene probe ...\n")
      
      isGeneProbe<-!(nm %in% grep("^A", nm, value=TRUE))
      geneProbeData<-subset(norm_all,isGeneProbe)
      rm(norm_all)
      
      cat("find largest IQR probe for each gene ...\n")
      
      arrayIQR<-apply(geneProbeData,1,IQR)
      probeName<-rownames(geneProbeData)
      largestProbeName<-findLargest(as.vector(probeName),arrayIQR,"hgu133a")
      
      probefile=paste0(qnormfile,".bnorm.largestprobe.RData")
      save(largestProbeName, file=probefile)
      
      cat("present each gene by largest IQR probe ...\n")
      
      nm<-rownames(geneProbeData)
      isLargestProbe<-nm %in% largestProbeName
      geneData<-subset(geneProbeData,isLargestProbe)
      rm(geneProbeData)
      geneNames<-unlist(mget(rownames(geneData),envir=hgu133aSYMBOL))
      rownames(geneData)<-geneNames
      
      save(geneData,file=paste0(bnormfile,".RData"))
      
      cat("calculate stdev of each gene ...\n")
      
      geneDataSd <- apply(geneData, 1, function(x){sd(x)})
      write.csv(x=geneDataSd, file=paste0(bnormfile, ".sd.csv"))
      
      cat("keep gene with larger stdev ...\n")
      
      kdata<-subset(geneData, geneDataSd > 0.8)
      rm(geneData)
      save(kdata, file=sd_rdata)
      
      cat(paste0("save result to ", sdfile, ".csv ...\n"))
      write.csv(x=kdata, file=paste0(sdfile, ".csv"))
    }else{
      cat(paste0("load data from ", sd_rdata, " ...\n"))
      load(sd_rdata)
    }
    
    cat("calculate optimal k-mean cluster count ...\n")
    kdata<-t(kdata)
    res <- kmeansGap(kdata, hclust=TRUE)
    save(res, file=kmeans_rdata)
  }
  
  cat("save image ...\n")
  png(file=paste0(sdfile, ".kmeans.gap.png"), width=4000, height=3000, res=300);
  plot.clusterGap(res)
  dev.off()
  
  cat("done\n")
}
