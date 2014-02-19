source("e:/sqh/programs/r/mirna/chris_v2/common.r")

# #read mirna data
# mirnaCount<-loadOriginalData("mirnaseq", "count")
# brcaIndex<-grep("brca_", colnames(mirnaCount), value=FALSE)
# brcaCount<-mirnaCount[,brcaIndex]
# saveData(brcaCount, "TCGA_mirnaseq_20130521_brca_count.csv")
mirnaCount<-loadData("TCGA_mirnaseq_20130521_brca_count.csv")

#keep only observed in at least 90% samples
samplecount<-ncol(mirnaCount) * 0.9
keep<-rowSums(mirnaCount > 0) >= samplecount
mirnaCount<-mirnaCount[keep,]

#rank
mirnaRank<-apply(mirnaCount, 2, function(x) rank(x,ties.method="average"))
targetMirnaRank<-mirnaRank[row.names(mirnaRank) %in% targetnames,]
targetMirnaRankOfSample<-apply(targetMirnaRank, 1, function(x) rank(x,ties.method="average"))

# #read rnaseq data
# rnaseqCount<-loadOriginalData("rnaseqv2", "count")
# brcaIndex<-grep("brca_", colnames(rnaseqCount), value=FALSE)
# brcaCount<-rnaseqCount[,brcaIndex]
# saveData(brcaCount, "TCGA_rnaseqv2_20130521_brca_count.csv")
rnaseqCount<-loadData("TCGA_rnaseqv2_20130521_brca_count.csv")

#keep only observed in at least 90% samples
samplecount<-ncol(rnaseqCount) * 0.9
keep<-rowSums(rnaseqCount > 0) >= samplecount
rnaseqCount<-rnaseqCount[keep,]

rnaseqRank<-apply(rnaseqCount, 2, function(x) rank(x,ties.method="average"))
rnaseqRankOfSample<-apply(rnaseqRank, 1, function(x) rank(x,ties.method="average"))

brca<-cbind(targetMirnaRankOfSample, rnaseqRankOfSample)
saveData(brca, "TCGA_mirnaseq_rnaseqv2_20130521_brca_rank_sample.csv")

