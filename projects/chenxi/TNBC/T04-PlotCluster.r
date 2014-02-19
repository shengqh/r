library("plotrix")
library("cluster")

plot.clusterGap <- function(x, ...){
  plot(1:nrow(x$gapStat), x$gapStat[,2], ylab = "gap statistic", xlab = "number of clusters", type = "l")
  plotrix::plotCI(1:nrow(x$gapStat), x$gapStat[,2], x$gapStat[,3], add = TRUE)
  title("Gap statistic for clustering")
}

setwd("/scratch/cqs/shengq1/breastcancer/final/trainingdataset_old_repeat")

load("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.RData");

png(file="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.gap.png", width=4000, height=3000, res=300);
plot.clusterGap(res)
dev.off()

daisyfile <- "breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.daisy.RData"

if (file.exists(daisyfile)){
  load(daisyfile);
}else{
  load("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.RData")

  kdata<-t(kdata)

  dissE <- daisy(kdata)

  save(dissE, file=daisyfile)
}

sk <- silhouette(res$cluster$cluster, dissE)

png(file="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.cluster.png", width=4000, height=3000, res=300);
plot(sk)
dev.off()
