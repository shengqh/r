setwd("E:/sqh/Dropbox/Projects/CNV/paper/")

library(ggplot2)
library(scales)
library(plyr)

data<-read.table("overlap_0.5.tsv",sep="\t",header=T)

df <- ddply(data, .(Method, Type), transform, 
            cum.perc = Reduce('+', list(Percentage/2,cumsum(c(0,head(Percentage,-1))))))

df$method2<-factor(df$Method, levels = c("cn.MOPS","CoNIFER","ExomeDepth","exomeCopy","Control-FREEC","CNVnator"))

jpeg(filename="figure3.jpg",width=2000, height=3000,res=300)

plot.new()

p <- ggplot(df, aes(x=method2, y=Percentage, fill=Class))
p <- p + geom_bar(stat = "identity", position="stack")
p <- p + xlab("") + ylab("Proportion")
p <- p + geom_text(aes(x=Method, y=cum.perc, ymax=Percentage, label=Count, hjust = 0.5, vjust=0.45))
p <- p + facet_wrap(~Type, ncol=1, scales="free_y")
p <- p + theme(axis.text.x=element_text(angle=90, vjust = 0.5)) 
p

grid.text(label = c("(A)","(B)","(C)"), x = 0.004, y = c(0.91, 0.63, 0.35), just = c(0,0), gp = gpar(col = "black", cex = 1.5, fontface = "plain", fontfamily = "serif"))

dev.off()

