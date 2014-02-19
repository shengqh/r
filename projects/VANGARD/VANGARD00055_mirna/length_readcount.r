len<-read.csv("H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/paper/bowtie2_hg19_mapping.csv")
png(filename="H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/paper/bowtie2_hg19_mapping.png", width=3000, height=2000, res=300)
boxplot(MatchedCount~Length, data=len)
dev.off()


len<-read.table("H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/paper/2516-01.bam.count.matched.len", header=T)
png(filename="H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/paper/2516-01.bam.count.matched.png", width=3000, height=2000, res=300)
boxplot(MatchedCount~Length, data=len)
dev.off()

len<-read.table("E:/sqh/Dropbox/papers/2013-QC-review/lan_readcount.txt", header=T)
png(filename="E:/sqh/Dropbox/papers/2013-QC-review/lan_readcount.txt.png", width=3000, height=2000, res=300)
len["logTotalReads"]<-log10(len$TotalReads)
boxplot(logTotalReads~LaneID, data=len, xlab="Lane", ylab="log10(Total reads)")
dev.off()
