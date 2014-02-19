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

drawpca<-function(result, res){
  clst<-res$cluster$cluster
  ncluster<-nrow(res$cluster$centers);
  
  colors=rainbow(ncluster)
  pointnames=c(97:(97+ncluster-1))
  aa<-rawToChar(as.raw(pointnames))
  clunames<-substring(aa,seq(1,nchar(aa)),seq(1,nchar(aa)))
  
  for (i in 1:ncluster){
    clunames[i] <- paste0(clunames[i], " (", length(clst[clst==i]), ")");
  }
  
  plot(result$x[,1],result$x[,2],xlab="PC1",ylab="PC2",main="PCA", col=colors[clst],type="n")
  points(result$x[,1],result$x[,2],col=colors[clst],pch=pointnames[clst])
  legend("topright",ncluster,legend=clunames,lty=1, col=colors)
}

iswin = (Sys.info()['sysname'] == "Windows");

if (iswin){
  root <-"D:/projects/BreastCancer/trainingdataset_rdata";
}else{
  root <-"/scratch/cqs/shengq1/breastcancer/final/trainingdataset_rdata";
}

setwd(root)
#original data
datafile<-"breastcancer_affymetrix.tsv";

function<-build_dataset(datafile, istraining=T, largestprobefile=""){
#quantile normalization data
qnormfile<-paste0(datafile,".qnorm");
qnormdatafile<-paste0(qnormfile,".RData");
if(!file.exists(qnormdatafile)){
  #read/transfer data
  data<-read.table(datafile,header=T,row.name=1);
  tdata<-t(data);
  rm(data);

  #quantile normalization
  normalize.quantiles.robust(x=tdata,copy=F);
  qnormdata<-log2(tdata);
  rm(tdata);
          
  save(qnormdata,file=qnormdatafile);
  #write.csv(qnormdata,file=qnormfile);
}

#batch normalization data
bnormfile<-paste0(qnormfile,".bnorm");
bnormdatafile<-paste0(bnormfile,".RData");
if(!file.exists(bnormdatafile)){
  load(qnormdatafile);

  #batch normalization
  batchdefinitionfile<-paste0(datafile,".batchdefinition");
  batch<-read.table(batchdefinitionfile)
  foo<-lm(t(qnormdata) ~ as.factor(batch$V1))
  rm(qnormdata)
  bnormdata<-t(foo$res)
  rm(foo)

  save(bnormdata,file=bnormdatafile);
  #write.csv(bnormdata,file=bnormfile);
}

#from probe to genes
genesfile=paste0(bnormfile,".genes")
genesdatafile<-paste0(genesfile,".RData");
if(!file.exists(genesdatafile)){
  load(bnormdatafile);
  
  nm <- as.vector(sapply(rownames(bnormdata), function(x) substr(x,3,nchar(x)),simplify=TRUE))
  rownames(bnormdata)<-nm
  
  #filter gene probe
  isGeneProbe<-!(nm %in% grep("^A", nm, value=TRUE))
  geneProbeData<-subset(bnormdata,isGeneProbe)
  rm(bnormdata)
  
  if(istraining){
    #find largest IQR probe for each gene
    arrayIQR<-apply(geneProbeData,1,IQR)
    probeName<-rownames(geneProbeData)
    largestProbeName<-findLargest(as.vector(probeName),arrayIQR,"hgu133a")
    
    #save largestProbeName
    probefile=paste0(bnormfile,".largestprobe");
    save(largestProbeName, file=paste0(probefile,".RData"));
    #write.csv(largestProbeName, probefile);
  }
  else{
    load(largestprobefile);
  }
  
  #present each gene by largest IQR probe
  nm<-rownames(geneProbeData)
  isLargestProbe<-nm %in% largestProbeName
  geneData<-subset(geneProbeData,isLargestProbe)
  rm(geneProbeData)
  
  geneNames<-unlist(mget(rownames(geneData),envir=hgu133aSYMBOL))
  rownames(geneData)<-geneNames
  
  #save gene data
  save(geneData,file=genesdatafile);
  #write.csv(geneData, file=genesfile);
}

#calculate stdev of each gene
allsddatafile<-paste0(genesfile, ".sd.RData");
if(!file.exists(allsddatafile)){
  load(genesdatafile);
  geneDataSd <- apply(geneData, 1, function(x){sd(x)});
  save(geneDataSd, file=allsddatafile);
}

#keep only genes whose std >= 0.8
sdfile=paste0(genesfile,".sd0.8");
sddatafile=paste0(sdfile,".RData");
if(!file.exists(sddatafile)){
  load(genesdatafile);
  load(allsddatafile);
  
  #keep gene with larger stdev
  kdata<-subset(geneData, geneDataSd > 0.8);
  rm(geneData);
  
  #save
  save(kdata, file=sddatafile);
  #write.csv(kdata, file=sdfile);
}

#caluclate kmeans count
kmeansfile<-paste0(sdfile,".kmeans");
kmeansdatafile<-paste0(kmeansfile,".RData");
if(!file.exists(kmeansdatafile)){
  load(sddatafile);

  #calculate optimal k-mean cluster count
  kdata<-t(kdata)
  res <- kmeansGap(kdata, hclust=TRUE)
  save(res, file=kmeansdatafile)
  #write.csv(res, file=kmeansfile);
} 

#load kmeans data
load(kmeansdatafile);
clst<-res$cluster$cluster
ncluster<-nrow(res$cluster$centers);

#save image
kmeanspngfile<-paste0(kmeansfile, ".gap.png");
if(!file.exists(kmeanspngfile)){
  png(file=kmeanspngfile, width=4000, height=3000, res=300);
  plot.clusterGap(res)
  dev.off()
}

#get updown genelist and gene data
updown<-paste0(kmeansfile, ".updown");
updowngenelistfile<-paste0(updown,".genelist.RData");
updowngenesfile<-paste0(updown,".genes.RData");
if(!file.exists(updowngenelistfile)){
  load(genesdatafile);
  genecount = nrow(geneData)
  genenames <- rownames(geneData)
  samplenames <- colnames(geneData)

  PERCENTAGE = 0.2
  MINCOVERAGE = 0.7

  genes<-list()
  genelist<-list()
  
  for (i in 1:ncluster) {
    pid<-clst[clst==i]
    pid<-as.character(names(pid))
    sub<-geneData[,pid]
    #rank genes in each sample
    pen<-apply(sub, 2, rank)
  
    #normalize to 0~1
    pen<-pen / genecount
  
    #if the gene is in top 20%
    up <- (pen >= (1 - PERCENTAGE))
  
    #if the gene is in bottom 20%
    down <- (pen <= PERCENTAGE)
  
    aa<-list()
    for (j in 1:genecount) {
      #get gene ranks in current cluster
      genej<-up[j,]
    
      #get top 20% genes in the cluster
      upj<-genej[genej==TRUE]
      sigu = length(upj) / length(genej)
      if(sigu > MINCOVERAGE){
        genes[genenames[j]] = 1
        aa[genenames[j]] <- "up"
        next
      }
    
      #get bottom 20% genes in the cluster
      genej<-down[j,] 
      downj<-genej[genej==TRUE]
      sigd = length(downj) / length(genej)
      if(sigd > MINCOVERAGE){
        genes[genenames[j]] = 1
        aa[genenames[j]] <- "down"
      }
    }
    genelist[[i]] <- aa
  }

  #save gene list
  save(genelist, file=updowngenelistfile);

  #save gene expression data
  isUpDown<-genenames %in% names(genes)
  updownData<-subset(geneData, isUpDown)
  save(updownData, file=updowngenesfile);
}

#updown gene pca analysis
updownpcafile<-paste0(updown,".pca.RData");
if(!file.exists(updownpcafile)){
  load(updowngenesfile);
  result<-prcomp(t(updownData))
  save(result, file=updownpcafile);
}

#updown gene pca image
updownpcapngfile<-paste0(updown,".pca.png");
if(!file.exists(updownpcapngfile)){
  load(updownpcafile);

  png(filename=updownpcapngfile,width=4000, height=4000,res=300);
  drawpca(result, res);
  dev.off();
}

#whole gene pca analysis
genepcadatafile<-paste0(genesfile, ".pca.RData");
if(!file.exists(genepcadatafile)){
  load(genesdatafile);
  result<-prcomp(t(geneData))
  save(result, file=genepcadatafile);
}

#whole gene pca image
genepcapngfile<-paste0(genesfile, ".pca.png");
if(!file.exists(genepcapngfile)){
  load(genepcadatafile);

  png(filename=genepcapngfile,width=4000, height=4000,res=300);
  drawpca(result, res);
  dev.off();
}

