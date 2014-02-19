source("e:/sqh/programs/r/mirna/chris/common.r")

filteredfile<-"NCI60_microarray_CEL_filtered.tsv"
nci60filtered<-read.table(filteredfile, sep="\t", header=T, check.names=FALSE, row.names=1)

mirna90<-read.table("NCI60_miRNA_GSE26375_TCGA90.tsv", sep="\t", header=T, check.names=FALSE, row.names=1)

targetMirna90RankOrdered<-getTargetMirnaRankOrdered(mirna90)
saveData(targetMirna90RankOrdered, "NCI60_miRNA_GSE26375_TCGA90.target.rank.ordered.csv")

mirnacolors <- colorRampPalette(c("green", "black", "red"))(ncol(targetMirna90RankOrdered))
names(mirnacolors)<-colnames(targetMirna90RankOrdered)

common<-colnames(nci60filtered) %in% colnames(mirna90)
only<-colnames(nci60filtered)[!common]
print(paste0("No mirna data of ", only))

y<-as.matrix(nci60filtered[,common])

hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="average"); 
hc <- hclust(as.dist(1-cor(y, method="pearson")), method="average") 

mircols<-mirnacolors[colnames(y)]
clab  <- matrix(c(rep("white", ncol(y)), mircols), ncol=2, byrow=FALSE)
colnames(clab)<-c("", "miRNA")

png(file=paste0(filteredfile, ".png"), width=3000, height=3000,res=300)
heatmap.plus(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=hmcols, ColSideColors=clab, labRow=rep("", nrow(y))) 
dev.off()

