library(ggplot2)

setwd("E:/sqh/Dropbox/career/papers/draft/2014-iTRAQ")

data<-read.table("Response_1_ScoreImprovement.tsv", header=T,row.names=1)

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

postscript(file="Figure1_ScoreImprovement.eps", width=800, height = 800)

ma_xloc <-19
ma_yloc <-120
ma.cor <- data.frame(PlexType=c("iTRAQ4", "iTRAQ8"), Text=c("A","B"))
p1<-ggplot(data, aes(AverageScore, ScoreImprovement)) + 
  geom_point(shape = 1) + 
  facet_wrap(~ PlexType) + 
  geom_hline(yintercept = 0) +
  xlab("(Processed score + Unprocessed score) / 2") +
  ylab("Processed score - Unprocessed score") + 
  geom_text(data=ma.cor, aes(x=ma_xloc, y=ma_yloc, label=Text), 
            colour="black", parse=FALSE, size=8)

hist_xloc=-25
hist_yloc=0.042
hist.cor <- data.frame(PlexType=c("iTRAQ4", "iTRAQ8"), Text=c("C","D"))
p2<-ggplot(data, aes(x=ScoreImprovement)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=2,
                 colour="black", fill="white") +
  facet_wrap(~ PlexType) +
  xlab("Processed score - Unprocessed score") +
  ylab("Density") +
  geom_text(data=hist.cor, aes(x=hist_xloc, y=hist_yloc, label=Text), 
          colour="black", parse=FALSE, size=8)

multiplot(p1, p2)
          
dev.off()
