setwd("E:/sqh/Dropbox/Projects/CNV/paper/")

readtypedata<-function(filename, method){
  data<-read.table(filename,sep="\t", header=T)
  tt<-table(data$type)
  mt<-melt(tt,id.vars=1:2)
  mtt<-cbind(c(rep(method,nrow(mt))), mt)
  return(mtt)
}

readlendata<-function(filename, method){
  data<-read.table(filename,sep="\t", header=T)
  mtt<-cbind(c(rep(method,nrow(data))), data[,1:6])
  colnames(mtt)[1]<-"Method"
  return(mtt)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(ggplot2)
library(reshape2)
library(plyr)

names<-c("aCGH","cn.MOPS","CNVnator","CoNIFER","Control-FREEC", "exomeCopy", "ExomeDepth")
filenames<-c("cgh.txt","cnmops.txt","cnvnator.txt","conifer.txt","freec.txt","exomeCopy.txt","ExomeDepth.txt")

datatype<-lapply(c(1:7), FUN=function(x){
  return (readtypedata(filenames[x],names[x]))
})
typedata<-rbind(datatype[[3]], 
                datatype[[1]], 
                datatype[[5]], 
                datatype[[6]],
                datatype[[7]],
                datatype[[2]], 
                datatype[[4]]
)
colnames(typedata)<-c("Method","CNV","Count")

datalen<-lapply(c(1:7), FUN=function(x){
  return (readlendata(filenames[x],names[x]))
})
lendata<-rbind(datalen[[5]][,c(1,6)], 
               datalen[[3]][,c(1,6)], 
               datalen[[6]][,c(1,6)],
               datalen[[7]][,c(1,6)],
               datalen[[1]][,c(1,6)], 
               datalen[[4]][,c(1,6)], 
               datalen[[2]][,c(1,6)]
)
lendata[,2]<-log2(lendata[,2])
names(lendata)[2]<-"Length"

jpeg(filename="figure1.jpg",width=4000, height=2000,res=300)

df <- ddply(typedata, .(Method), transform, 
            cum.perc = Reduce('+', list(Count/2,cumsum(c(0,head(Count,-1))))))
df[13,4]<--50
df[14,4]<-210
g1<-ggplot(df, aes(x=Method, y=Count, fill=CNV))
g1<-g1+geom_bar(stat="identity", position="stack") + theme_bw() + xlab("")
g1 <- g1 + geom_text(aes(x=Method, y=cum.perc, ymax=cum.perc, label=Count, hjust = 0.5, vjust=0.45))
g1<-g1 +  theme(axis.text.x=element_text(angle=90, vjust = 0.5))

g2<-qplot(Method, Length, data=lendata, geom=c("boxplot"), fill=Method, ylab="log2(CNV length)", xlab="") + theme_bw() + 
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) 

multiplot(g1,g2,cols=2)

#text(x=bplt, y=lenarray+180, labels=as.character(lenarray))

grid.text(label = "(A)", x = 0.01, y = 0.95, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))
grid.text(label = "(B)", x = 0.5, y = 0.94, just = c(0,0), gp = gpar(col = "black", cex = 2, fontface = "plain", fontfamily = "serif"))

dev.off()

summary(2 ** (lendata[lendata$Method=="aCGH",2]))
