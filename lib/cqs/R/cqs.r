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

readCDT <- function(filename) {
  fname <- sub('.cdt$', '', filename) # get rid of the extension
  cdt <- read.table(filename, sep='\t', header=TRUE, row.names=NULL)
  
  # we only need the first column of the CDT file, which contains
  # the order information, and the third column, which contains the
  # labels. The rest of the file contains the data matrix.
  rown <- as.character(cdt[,"gene"])
  coln <- colnames(cdt)
  firstRow <- 1 + which(rown=="EWEIGHT")
  firstCol <- 1 + which(coln=="GWEIGHT")
  
  gid <- as.character(cdt[,"gene"])[firstRow:nrow(cdt)]
  aid <- cdt[rown=="AID",][firstCol:ncol(cdt)]
  aid <- as.character(as.matrix(aid))
  
  # names all start with 'GENE' or 'NODE' (or 'ARRY') and end with 'X'
  gene.order <- 1 + as.numeric(substring(gid, 5, nchar(gid)-1))
  arry.order <- 1 + as.numeric(substring(aid, 5, nchar(aid)-1))
  
  # Because Cluster reorders things and because hclust and plclust wants
  # to do the same, we have to reinvert the ordering during passage from
  # one to the other
  gene.labels <- as.character(cdt$NAME)[firstRow:nrow(cdt)][order(gene.order)]
  arry.labels <- coln[firstCol:ncol(cdt)][order(arry.order)]
  
  temp <- as.matrix(cdt[firstRow:nrow(cdt), firstCol:ncol(cdt)])
  temp <- temp[order(gene.order), order(arry.order)]
  data <- matrix(as.numeric(temp), ncol=ncol(temp))
  dimnames(data) <- list(gene.labels, arry.labels)
  
  # Columns 2 and 3 describe the two branches below each node.
  # Nodes are listed from bottom to top since clustering is
  # agglomerative.
  #
  foo <- function(alt) {
    # Again, we get the numeric part of the label
    base <- as.numeric(substring(alt, 5, nchar(alt)-1))
    # We also need to know whether the label is a "GENE" or a "NODE".
    # The 'hclust' objects use negative integers to indicate nodes.
    type1 <- rep(1, length(base))
    type1[substring(alt, 1, 4) %in% c('GENE', "ARRY")] <- -1
    base <- base*type1  # make nodes negative
    adder <- (type1-1)/2  # offset the negatives to change from starting
    # at 1 to starting at 0.
    base + adder
  }
  
  readdata <-function(datafilename, aOrder, aLabels){
    # The gtr file contains the "distances" in column 4. Actually,
    # Eisen's Cluster program reports similarities instead of
    # distances. This fix assumes that some kind of correlation was
    # the meaure of similarity....
    xtr <- read.table(datafilename, sep='\t',
                      header=FALSE, as.is=TRUE)
    d.height <- 1 - xtr$V4
    d.merge1 <- foo(xtr$V3)
    d.merge2 <- foo(xtr$V2)
    result <- list(merge=as.matrix(cbind(d.merge1, d.merge2)),
                   height=d.height,
                   order=aOrder,
                   labels=aLabels,
                   method='modified centroid',
                   call=NULL,
                   dist.method='Pearson correlation')
    class(result) <- 'hclust'
    result
  }
  
  gtrfile <-paste(fname, 'gtr', sep='.')
  if(file.exists(gtrfile))  {
    genedata <- readdata(gtrfile, gene.order, gene.labels);
  }
  else{
    genedata <- 0;
  }
  
  atrfile <-paste(fname, 'atr', sep='.')
  if (file.exists(atrfile))    {
    arrydata <- readdata(atrfile, arry.order, arry.labels);
  }
  else{
    arrydata <- 0;
  }
  
  list(gene=genedata,arry=arrydata, data=data)
}
