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
  
  temp1 <- log(sum(by(data, factor(class), pw.dist)))q
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

drawpca<-function(pca_res, kmeans_res){
  clst<-kmeans_res$cluster$cluster
  ncluster<-nrow(kmeans_res$cluster$centers);
  
  colors=rainbow(ncluster)
  pointnames=c(97:(97+ncluster-1))
  aa<-rawToChar(as.raw(pointnames))
  clunames<-substring(aa,seq(1,nchar(aa)),seq(1,nchar(aa)))
  
  for (i in 1:ncluster){
    clunames[i] <- paste0(clunames[i], " (", length(clst[clst==i]), ")");
  }
  
  plot(pca_res$x[,1],pca_res$x[,2],xlab="PC1",ylab="PC2",main="PCA", col=colors[clst],type="n")
  points(pca_res$x[,1],pca_res$x[,2],col=colors[clst],pch=pointnames[clst])
  legend("topright",ncluster,legend=clunames,lty=1, col=colors)
}

#quantile normalization data
quantile_normalization<-function(myfile){
  myfile$qnormfile<-paste0(myfile$datafile,".qnorm");
  myfile$qnormdatafile<-paste0(myfile$qnormfile,".RData");
  if(!file.exists(myfile$qnormdatafile)){
    #read/transfer data
    data<-read.table(myfile$datafile,header=T,row.name=1);
    tdata<-t(data);
    rm(data);
    
    #quantile normalization
    normalize.quantiles.robust(x=tdata,copy=F);
    qnormdata<-log2(tdata);
    rm(tdata);
    
    save(qnormdata,file=myfile$qnormdatafile);
    #write.csv(qnormdata,file=qnormfile);
    rm(qnormdata);
  }
  return(myfile);
}

#batch normalization data
batch_normalization<-function(myfile){
  myfile$bnormfile<-paste0(myfile$qnormfile,".bnorm");
  myfile$bnormdatafile<-paste0(myfile$bnormfile,".RData");
  if(!file.exists(myfile$bnormdatafile)){
    load(myfile$qnormdatafile);
    
    #batch normalization
    batchdefinitionfile<-paste0(myfile$datafile,".batchdefinition");
    batch<-read.table(batchdefinitionfile)
    foo<-lm(t(qnormdata) ~ as.factor(batch$V1))
    rm(qnormdata)
    bnormdata<-t(foo$res)
    rm(foo)
    
    save(bnormdata,file=myfile$bnormdatafile);
    #write.csv(bnormdata,file=bnormfile);
    rm(bnormdata);
  }
  return (myfile);
}

#from probe to genes
probe_to_genes<-function(myfile, largeprobefile=""){
  myfile$genesfile=paste0(myfile$bnormfile,".genes")
  myfile$genesdatafile<-paste0(myfile$genesfile,".RData");

  if(!file.exists(myfile$genesdatafile)){
    load(myfile$bnormdatafile);
  
    nm <- as.vector(sapply(rownames(bnormdata), function(x) substr(x,3,nchar(x)),simplify=TRUE))
    rownames(bnormdata)<-nm
  
    #filter gene probe
    isGeneProbe<-!(nm %in% grep("^A", nm, value=TRUE))
    geneProbeData<-subset(bnormdata,isGeneProbe)
    rm(bnormdata)
  
    #get largestProbeName
    if(largeprobefile == ""){
      #find largest IQR probe for each gene
      arrayIQR<-apply(geneProbeData,1,IQR)
      probeName<-rownames(geneProbeData)
      largestProbeName<-findLargest(as.vector(probeName),arrayIQR,"hgu133a")
  
      #save largestProbeName
      probefile=paste0(myfile$bnormfile,".largestprobe");
      myfile$probedatafile=paste0(probefile,".RData");
      save(largestProbeName, file=myfile$probedatafile);
      #write.csv(largestProbeName, probefile);
    }
    else{
      load(largeprobefile);
    }
  
    #present each gene by largest IQR probe
    nm<-rownames(geneProbeData)
    isLargestProbe<-nm %in% largestProbeName
    rm(largestProbeName)
    
    geneData<-subset(geneProbeData,isLargestProbe)
    rm(geneProbeData)
  
    geneNames<-unlist(mget(rownames(geneData),envir=hgu133aSYMBOL))
    rownames(geneData)<-geneNames
  
    #save gene data
    save(geneData,file=genesdatafile);
    #write.csv(geneData, file=genesfile);
    
    rm(geneData);
  }
  return (myfile);
}

#calculate stdev of each gene
calculate_gene_stdev<-function(myfile){
  myfile$allsddatafile<-paste0(myfile$genesfile, ".sd.RData");
  
  if(!file.exists(myfile$allsddatafile)){
    load(myfile$genesdatafile);
    geneDataSd <- apply(geneData, 1, function(x){sd(x)});
    save(geneDataSd, file=myfile$allsddatafile);
    rm(geneDataSd);
  }
  
  return (myfile);
}

#keep only genes whose std >= 0.8
filter_gene_by_stdev<-function(myfile){
  myfile$sdfile=paste0(myfile$genesfile,".sd0.8");
  myfile$sddatafile=paste0(myfile$sdfile,".RData");
  if(!file.exists(myfile$sddatafile)){
    load(myfile$genesdatafile);
    load(myfile$allsddatafile);
    
    #keep gene with larger stdev
    kdata<-subset(geneData, geneDataSd > 0.8);
    rm(geneData);
    
    #save
    save(kdata, file=myfile$sddatafile);
    #write.csv(kdata, file=sdfile);
    
    rm(kdata)
  }
  
  return (myfile);
}

#caluclate kmeans
calculate_kmeans(myfile){
  myfile$kmeansfile<-paste0(myfile$sdfile,".kmeans");
  myfile$kmeansdatafile<-paste0(myfile$kmeansfile,".RData");
  if(!file.exists(myfile$kmeansdatafile)){
    load(myfile$sddatafile);
    
    #calculate optimal k-mean cluster count
    kdata<-t(kdata)
    res <- kmeansGap(kdata, hclust=TRUE)
    save(res, file=myfile$kmeansdatafile)
    #write.csv(res, file=kmeansfile);
    
    rm(res);
  } 
  return (myfile);
}

#save kmeans image
save_kmeans_image(myfile){
  myfile$kmeanspngfile<-paste0(myfile$kmeansfile, ".gap.png");
  if(!file.exists(myfile$kmeanspngfile)){
    load(myfile$kmeansdatafile);
    png(file=myfile$kmeanspngfile, width=4000, height=3000, res=300);
    plot.clusterGap(res)
    dev.off()
    rm(res)
  }
  
  return (myfile);
}

#get updown genelist and gene data
get_updown_data(myfile){
  myfile$updown<-paste0(myfile$kmeansfile, ".updown");
  myfile$updowngenelistfile<-paste0(myfile$updown,".genelist.RData");
  myfile$updowngenesfile<-paste0(myfile$updown,".genes.RData");
  if(!file.exists(myfile$updowngenelistfile)){
    load(myfile$genesdatafile);
    load(myfile$kmeansdatafile);
    
    
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
}



iswin = (Sys.info()['sysname'] == "Windows");

if (iswin){
  root <-"D:/projects/BreastCancer/trainingdataset_rdata";
}else{
  root <-"/scratch/cqs/shengq1/breastcancer/final/trainingdataset_rdata";
}

setwd(root)

#original data
myfile<-list();
myfile$datafile<-"breastcancer_affymetrix.tsv";

myfile<-quantile_normalization(myfile);
myfile<-batch_normalization(myfile);
myfile<-probe_to_genes(myfile);
myfile<-calculate_gene_stdev(myfile);
myfile<-filter_gene_by_stdev(myfile);
myfile<-calculate_kmeans(myfile);
myfile<-save_kmeans_image(myfile);


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

