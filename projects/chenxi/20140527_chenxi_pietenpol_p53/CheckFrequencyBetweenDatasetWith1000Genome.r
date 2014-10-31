setwd("H:/shengquanhu/projects/chenxi/20140527_chenxi_pietenpol_p53/impute2")
data<-read.table("Pietenpol_p53.genotype.info", header=T, row.names=1, sep="\t")

subdata<-data[data$X1000GenomeAllele1 != " ",]
plot(subdata$Allle2Frequency, subdata$X1000GenomeAllele2Frequency)
