library("preprocessCore")

rootdir<-"E:/sqh/Dropbox/papers/2013-RNAseqMicroarrayComparison";

datafile<-paste0(rootdir, "/Figure3_BRCA.tsv");
data<-read.table(datafile, header=T, row.names=1, check.names=F, sep="\t")

#normalize microarray data by quantiles
microdata<-data[,1:106]
microdata[is.na(microdata)]<-1.0
dd<-as.matrix(microdata)
ddd<-normalize.quantiles.robust(x=dd,copy=T)

#paired t-test
data$AgilentFold = apply(ddd, 1, function(x) log2(median(exp(x[54:106])) / median(exp(x[1:53])) ))
data$AgilentPValue = apply(ddd, 1, function(x) t.test(x[1:53], x[54:106], paired=T)$p.value)
data$AgilentFDR = p.adjust(data$AgilentPValue, method="fdr")
res<-data[,(ncol(data)-2):ncol(data)];
sdata<-res[order(rownames(res)),]

#read rnaseq test result from count data
rpkm<-read.csv(paste0(rootdir, "/deseq.csv"), row.names=1)[,c(1,6,7,8)]
rownames(rpkm)<-rpkm[,1]
rpkm<-rpkm[,c(2,3,4)]
colnames(rpkm)<-c("RPKMFold","RPKMPValue","RPKMFDR")
srpkm<-rpkm[order(rownames(rpkm)),]

#read raw data
rdata<-read.table(paste0(rootdir, "/Figure3_BRCA.tsv"), header=T, row.names=1, check.names=F, sep="\t")
rsdata<-rdata[order(rownames(rdata)),]
rsdata$RPKMLog2Mean<-apply(rsdata,1,function(x) log2(mean(x[213:318])))

#output fdr result
mdata<-cbind(sdata, srpkm, rsdata$RPKMLog2Mean)
colnames(mdata)[ncol(mdata)]<-"RPKMLog2Mean"

write.csv(mdata, file=paste0(rootdir, "/Figure3_BRCA_PValue.csv"))


