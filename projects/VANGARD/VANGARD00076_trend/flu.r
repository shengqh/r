library("Kendall")

setwd("D:/projects/Flu")

fpkmfiles <- list.files(pattern="*.fpkm", full.names=TRUE)

trendfunc <- function(x){
  if (x[3] == 0 && x[4] == 0 && x[5] == 0 && x[6] == 0){
    return (1)
  }
  else{
    v <- c(as.double(x[3]), as.double(x[4]),as.double(x[5]),as.double(x[6]))
    mv <-MannKendall(v)
    return (as.double(mv$sl))
  }
}

sdfunc <- function(x){
  if (x[3] == 0 && x[4] == 0 && x[5] == 0 && x[6] == 0){
    return (0)
  }
  else{
    v <- c(as.double(x[3]), as.double(x[4]),as.double(x[5]),as.double(x[6]))
    return (sd(v))
  }
}

rsdfunc <- function(x){
  if (x[3] == 0 && x[4] == 0 && x[5] == 0 && x[6] == 0){
    return (0)
  }
  else{
    v <- c(as.double(x[3]), as.double(x[4]),as.double(x[5]),as.double(x[6]))
    m <- mean(v)
    return (sd(v) * 100.0 / m)
  }
}

for (fpkmfile in fpkmfiles){
  fpkm <- read.table(fpkmfile, sep="\t")
  colnames(fpkm) = c("gene_short_name",  "locus",	"FPKM_0",	"FPKM_1","FPKM_3",	"FPKM_7")
  fpkm2 <- subset(fpkm, fpkm["FPKM_0"] != "FPKM")
  Pvalue<-apply(fpkm2,1,trendfunc)
  Stdev<-apply(fpkm2,1,sdfunc)
  Rsd<-apply(fpkm2,1,rsdfunc)
  fpkm3<-cbind(fpkm2, Pvalue, Stdev, Rsd)
  write.table(fpkm3, file=paste0(fpkmfile,".trend"), row.names=FALSE, sep="\t", quote=FALSE)
}