Goal: find where tcga laml RNASeq V1 data ranks SMC2 and SMC3 in its list of 
# genes differentially expressed in STAG2 low expressors vs high expressors. Use
# RPKM data in bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3.1.7.0.tar.gz found in 
# https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/laml/cgcc/bcgsc.ca/illuminaga_rnaseq/rnaseq/

# Note: to install bioconductor packages, use code like this
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
# instead of install.packages()

# 1. 30 points 

# Create a Biobase ExpressionSet object called eset1n (eset using RNASeq Version
# 1, normalized) out of *.gene.quantification.txt files (e.g. 
# TCGA-AB-2803-03A-01T-0734-13.gene.quantification.txt) in your 
# /data/laml/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3.1.7.0 folder. Use the RPKM 
# column as the normlized values. Do not include inv16, t8.21, bcr-abl (ba) or
# PML-RARA AML patients in eset1n (hint, it may be easier to take them out after
# they first get into the ExpressionSet).

###!!! set as your own directory and remove the following line

library("Biobase")

setwd("H:/shengquanhu/projects/TCGA/laml/data/rnaseq/illuminaga_rnaseq/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3.1.7.0")
files<-list.files(pattern="*.gene.quantification.txt")
patients<-substr(files, 1,12)

genes<-row.names(read.table(files[1], sep="\t", header=T, row.names=1, check.names=F))

data<-matrix(nrow=length(genes), ncol=length(patients), dimnames=list(genes, patients))
for(index in c(1:length(files))){
  data[,index]<-read.table(files[index], sep="\t", header=T, row.names=1)$RPKM
}

eset1n<-ExpressionSet(data)
###!!! I don't know how to remove inv16, t8.21, bcr-abl (ba) or PML-RARA AML patients in eset1n, is there any file indicating this information?

# 2. 30 points 

# Using eset1n from #1, plot a histogram of STAG2 values with a vertical line
# at a picked cut-off for STAG2 low expressors. Place the IDs of low expressors
# in a vector of strings called low1 and create a vector of strings called hi1
# of IDs of fourth quartile STAG2 expressors. Use low1 and hi1 to reduce eset1n
# to only patients in these two groups. Finally, make a data.frame with
# one column, stag2low, that contains one for patients low in Stag2 and zero
# otherwise, and add it to the eset1n object using pData() from the Biobase package.

genenames<-unlist(lapply(rownames(exprs(eset1n)), function(x){
  m<-regexpr("[^|]*", x, perl=TRUE)
  unlist(regmatches(x, m))
}))

stag2<-exprs(eset1n)[genenames=="STAG2",]
hist(stag2, breaks=20)
quan<-quantile(stag2)
abline(v=quan[2], col="red")

low1<-names(stag2)[stag2 < quan[2]]
hi1<-names(stag2)[stag2 > quan[4]]

low1data<-exprs(eset1n)[, names(stag2) %in% low1]
hi1data<-exprs(eset1n)[, names(stag2) %in% hi1]
lowhidata<-cbind(low1data, hi1data)
eset1n<-ExpressionSet(lowhidata)
stag2low<-data.frame(stag2low=c(rep(1, length(low1)), rep(0, length(hi1))))
pData(eset1n)<-stag2low

# 3. 40 points 

# Run voom() followed by lmFit() followed by eBayes() in the bioconductor 
# package limma on eset1n to find genes differentially expressed in low versus 
# high STAG2 expressors. With genes ordered by increasing BH adjusted P-values,
# where do SMC2 and SMC3 appear in list? Where did they appear in the list that
# stag2SigJoint.R produced using U133 joint AML-MDS data?

library(limma)

vdata<-voom(exprs(eset1n))
ldata<-lmFit(vdata)
edata<-eBayes(ldata)
pv<-data.frame(pvalue=edata$p.value, BH=p.adjust(edata$p.value))
pv<-pv[order(pv$BH),]

pvgenes<-unlist(lapply(rownames(pv), function(x){
  m<-regexpr("[^|]*", x, perl=TRUE)
  unlist(regmatches(x, m))
}

which(pvgenes %in% "SMC2")
which(pvgenes %in% "SMC3")

