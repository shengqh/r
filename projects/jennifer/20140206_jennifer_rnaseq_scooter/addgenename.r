setwd("H:/shengquanhu/projects/Jennifer/20140206_jennifer_rnaseq_scooter/")

map<-read.table("H:/shengquanhu/projects/database/ensembl_gtf/Homo_sapiens.GRCh37.74_chr1-22-X-Y-M.map", header=T, row.names=1)

data<-read.csv("MultiRankSeq2.html.diff.csv", header=T, row.names=1, check.names=F)
cdata<-cbind(map[rownames(data),], data)
colnames(cdata)[1]<-"Gene"
write.csv(cdata, "4v4.gene.diff.csv")

data<-read.csv("3v3.diff.count.csv", header=T, row.names=1, check.names=F)
cdata<-cbind(map[rownames(data),], data)
colnames(cdata)[1]<-"Gene"
write.csv(cdata, "3v3.gene.diff.csv")
