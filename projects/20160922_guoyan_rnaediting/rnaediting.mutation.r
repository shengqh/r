library(ggplot2)
setwd("H:/shengquanhu/projects")
dat=read.table("rnaediting.mutation.tsv", sep="\t", header=T)
numberOfSamples=length(unique(dat$File))
sqrtOfSample = ceiling(sqrt(numberOfSamples))

png(file="rnaediting.mutation.png", width=3000, height=3000, res=300)
g=ggplot(data=dat, aes(File, Percentage, fill=factor(Mutation))) + 
  geom_bar(stat="identity") + coord_flip()
print(g)
dev.off()
