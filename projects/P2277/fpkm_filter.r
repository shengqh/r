setwd("D:/projects/P2277/Cluster");

geneData <- read.table("fpkm.txt", header=T)

len<-length(colnames(geneData));

geneData[2:len] <- log(geneData[2:len] + 1);

logGeneDataSd <- apply(geneData[2:len], 1, function(x){sd(x)})

sortedsd<-sort(logGeneDataSd)

top=0.05;

topindex = trunc((1-top) * length(sortedsd))

topsd<-sortedsd[topindex]

kdata<-subset(geneData, logGeneDataSd > topsd)

write.table(x=kdata,  file="fpkm_0.05.txt", sep="\t", quote=FALSE, row.names=FALSE);
