setwd("H:/shengquanhu/projects/vickers/20150706_parclip_gsnap_3018-KCV-15/gsnap_smallRNA_count/result")
data<-read.table("t2c.tsv",sep="\t",header=T)

library(ggplot2)

png("t2c_rate_in_unique_reads.png", width=2000,height=2000,res=200)
ggplot(data, aes(Category, UniqueT2CRate, color=Category)) + geom_violin() + facet_grid(File ~ .) + ylab("T2C rate in unique reads")
dev.off()

png("t2c_rate_in_all_reads.png", width=2000,height=2000,res=200)
ggplot(data, aes(Category, TotalT2CRate, color=Category)) + geom_violin() + facet_grid(File ~ .) + ylab("T2C rate in all reads")
dev.off()

png("t2c_averageSitePer10Bases_in_unique_reads.png", width=2000,height=2000,res=200)
ggplot(data, aes(Category, AvergeT2CIn10BasesOfUniqueRead, color=Category)) + geom_violin() + facet_grid(File ~ .) + ylab("Average T2C mutation sites per 10 bases in unique reads")
dev.off()

png("t2c_averageSitePer10Bases_in_all_reads.png", width=2000,height=2000,res=200)
ggplot(data, aes(Category, AverageT2CIn10BasesOfTotalRead, color=Category)) + geom_violin() + facet_grid(File ~ .) + ylab("Average T2C mutation sites per 10 bases in all reads")
dev.off()
