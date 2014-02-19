setwd("/scratch/cqs/shengq1/breastcancer/final/trainingdataset_old_repeat")

load("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.RData")

load("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.RData")

clst<-res$cluster$cluster

#pca analysis
result<-prcomp(t(geneData))

#draw image
colors=rainbow(14)

pointnames=c(97:110)

aa<-rawToChar(as.raw(pointnames))

clunames<-substring(aa,seq(1,nchar(aa)),seq(1,nchar(aa)))

png(filename="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.pca.png",width=4000, height=4000,res=300);

plot(result$x[,1],result$x[,2],xlab="PC1",ylab="PC2",main="PCA", col=colors[clst],type="n")

points(result$x[,1],result$x[,2],col=colors[clst],pch=pointnames[clst])

legend("topright",14,legend=clunames,lty=1, col=colors)

dev.off();
