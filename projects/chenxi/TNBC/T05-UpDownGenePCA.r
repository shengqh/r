iswin = (Sys.info()['sysname'] == "Windows")

if (iswin)
  setwd("D:/projects/BreastCancer/final")

if(!iswin)
  setwd("/scratch/cqs/shengq1/breastcancer/final")

load("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.RData")

load("breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.RData")

clst<-res$cluster$cluster

genecount = nrow(geneData)

genenames <- rownames(geneData)

samplenames <- colnames(geneData)

PERCENTAGE = 0.2

MINCOVERAGE = 0.7

genes<-list()

genelist<-list()

for (i in 1:14) {
  #i=1
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
    #j=1
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

save(genelist, file="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.updown.genelist.RData");

#filter out updown gene data matrix
isUpDown<-genenames %in% names(genes)

updownData<-subset(geneData, isUpDown)
save(updownData, file="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.updown.genes.RData");

#pca analysis
result<-prcomp(t(updownData))
save(result, file="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.updown.pca.RData");

#draw image
colors=rainbow(14)

pointnames=c(97:110)

aa<-rawToChar(as.raw(pointnames))

clunames<-substring(aa,seq(1,nchar(aa)),seq(1,nchar(aa)))

png(filename="breastcancer_affymetrix.tsv.qnorm.bnorm.genes.sd0.8.kmeans.updown.pca.png",width=4000, height=4000,res=300);

plot(result$x[,1],result$x[,2],xlab="PC1",ylab="PC2",main="PCA", col=colors[clst],type="n")

points(result$x[,1],result$x[,2],col=colors[clst],pch=pointnames[clst])

legend("topright",14,legend=clunames,lty=1, col=colors)

dev.off();
