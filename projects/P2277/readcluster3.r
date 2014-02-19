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

if(1) {
    filename <- "D:/projects/P2277/Cluster/fpkm_0.05.cdt";
    if (!file.exists(filename)){
      next;
    }
    
    cdt <- readCDT(filename)
    a <- cdt$arry
    
    fname <- sub('.cdt$', '', filename) # get rid of the extension
    pngname<-paste0(fname,".png");
    
    png(filename=pngname, width=2000, height=3000,res=300);
    labels = rep("", length(a$labels))
    par(mar=c(10,5,2,1))
    plot(a, labels = labels, cex = 1, yaxt="n", ylab="P2277 top 0.05 sd(log(fpkm+1))", main="")
    mtext(a$labels, side = 1, line = -1, at = order(a$order), las = 2, out=FALSE)
    dev.off();  
}