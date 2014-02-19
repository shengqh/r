library ("affy")
library ("limma")

cqsMergeAffyBatch<-function (x, y, annotation = paste(annotation(x), annotation(y)), 
                             description = NULL, notes = character(0), ...) 
{
  adim <- dim(intensity(x))[1]
  if ((nrow(x) != nrow(y)) || (ncol(x) != ncol(y))) 
    stop("cannot merge chips of different sizes !")
  if (cdfName(x) != cdfName(y)) 
    warning("cdfName mismatch (using the cdfName of x)!")
  if (is.null(description)) {
    description <- new("MIAME")
    description@title <- "Merged AffyBatch"
  }
  lx <- length(x)
  ly <- length(y)
  phenodata <- phenoData(x)
  pData(phenodata) <- rbind(pData(x), pData(y))
  protocoldata <- protocolData(x)
  pData(protocoldata) <- rbind(pData(protocolData(x)), pData(protocolData(y)))
  notes(description) <- "Merged";

  #if merge two samples with same sampleName, sampleNames(phenodata) will be automatically
  #changed to different names, which will cause unequalment between sampleNames of data and 
  #sampleNames of phenoData.
  myexprs<-cbind(intensity(x), intensity(y));
  colnames(myexprs)<-sampleNames(phenodata);
  
  return(new("AffyBatch", exprs = myexprs, 
             phenoData = phenodata, experimentData = description, 
             cdfName = cdfName(x), nrow = nrow(x), ncol = ncol(x), 
             annotation = x@annotation, protocolData = protocoldata))
}

#covdesc<-"I:\\projects\\BreastCancer\\DifferentChipTypes\\covdesc.txt";
cqsReadAffyMixed<-function(covdesc="covdesc.txt") {
  cat("reading sample file names ...\n");
  
  samples <- read.AnnotatedDataFrame(covdesc, sep = "")
  
  files.to.read <- rownames(pData(samples))
  
  cat("\t", files.to.read, "\n", sep ="\n\t");
  
  cat("reading sample data ...\n");
  
  esets <- lapply(files.to.read, function(x) {
    cat("\treading sample : ", x, "\n", sep= " ");
    ReadAffy(filenames = x);
  })

  nms <- sapply(esets, cdfName)
  u.nms <- unique(nms)
  merged <- list()

  cat("\nmerging samples from same chip type ...\n");
  
  for (i in 1:length(nms)) {
    cat("\tprocessing",i,"\n", sep=" ");
    
    so.far <- merged[nms[i]]
    if (is.null(so.far[[1]])) {
      merged[[nms[i]]] <- esets[[i]]
    }
    else {
      merged[[nms[i]]] <- merge.AffyBatch(so.far[[1]], 
                                          esets[[i]])
    }
  }
  
  rm(esets);
  gc();

  nms <- unique(sapply(merged, cdfName))
  acc <- list(probeNames(merged[[1]]))
  nm <- cdfName(merged[[1]])[1]
  idxs <- list()
  idxs[nm] <- acc
  common <- as.vector(unlist(acc))
  for (i in 2:length(merged)) {
    nm <- cdfName(merged[[i]])[1]
    if (is.null(idxs[nm][[1]])) {
      acc <- list(probeNames(merged[[i]]))
      idxs[nm] <- acc
      common <- intersect(common, as.vector(unlist(acc)))
    }
  }
  smallest <- 1
  for (i in 2:length(merged)) {
    if (length(probeNames(merged[[i]])) < length(probeNames(merged[[smallest]]))) {
      smallest <- i
    }
  }
  all.smallest <- unique(probeNames(merged[[smallest]]))
  template <- merged[[smallest]][,1]
  not.shared <- setdiff(all.smallest, common)
  idxs <- indexProbes(template, "both")
  e <- exprs(template)
  e[unlist(idxs[is.element(all.smallest, not.shared)]), ] <- 0
  exprs(template) <- e

  cat("\nallocate memory for merged result ...\n");
  
  result <- template
  for (i in 2:length(files.to.read)) {
    result <- cqsMergeAffyBatch(result, template)
  }
  to.change <- exprs(result)
  idxs.small <- unlist(idxs[common])
  o <- order(names(idxs.small))
  idxs.small <- idxs.small[o]
  changing <- 1
  celnames<-rep('',length(files.to.read));

  cat("\ncopying data ...\n");
  
  for (i in 1:length(merged)) {
    idxs.large <- unlist(indexProbes(merged[[i]], "both")[common])
    o <- order(names(idxs.large))
    idxs.large <- idxs.large[o]
    for (j in 1:length(merged[[i]])) {
      to.copy <- exprs(merged[[i]][,j]);
      to.change[idxs.small, changing] <- to.copy[idxs.large, 1];
      celnames[changing] <- basename(sampleNames(merged[[i]][,j]));
      changing <- changing + 1;
    }
  }
  exprs(result) <- to.change
  
  cat("\nmerge finished.\n");
  
  return (result);
}

#result<-cqsReadAffyMixed(covdesc="I:\\projects\\BreastCancer\\Dataset\\celfile.list")
covdesc="I:\\projects\\BreastCancer\\DifferentChipTypes\\covdesc1.txt";
result<-cqsReadAffyMixed(covdesc=covdesc);
result.rma<-rma(result)
write.exprs(x=result.rma, file="I:\\projects\\BreastCancer\\Dataset\\breast_cancer_affy.tsv")
