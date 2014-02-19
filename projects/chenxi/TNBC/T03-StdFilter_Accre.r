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

## plot method for objects of class "clusterGap"
plot.clusterGap <- function(x, ...){
  plot(1:nrow(x$gapStat), x$gapStat[,2], ylab = "gap statistic", xlab = "number of clusters", type = "l")
  plotrix::plotCI(1:nrow(x$gapStat), x$gapStat[,2], x$gapStat[,3], add = TRUE)
  title("Gap statistic for clustering")
}

setwd("/scratch/cqs/shengq1/breastcancer/final/trainingdataset_old_repeat")

load("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.RData")

geneData.sd <- apply(geneData, 1, function(x){sd(x)})

write.csv(x=geneData.sd, file="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd.csv")

kdata<-subset(geneData, geneData.sd > 0.8)

rm(geneData)

save(kdata, file="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.RData")

#load("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.RData")

write.table(x=kdata, file="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.tsv",sep="\t")

kdata<-t(kdata)

res <- kmeansGap(kdata, hclust=TRUE)

save(res, file="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.RData")
