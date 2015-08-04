
setwd("H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/")  

comparisons=list(
  "HiSeq" = c("H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/hiseq/star_genetable/result/FFPE_FF_HiSeq_gene.count",
              "H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/hiseq/star_deseq2/result/HiSeq_FFPE_VS_FF.design", 
              "HiSeq_FF", "HiSeq_FFPE"),
  "MiSeq" = c("H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/miseq/star_genetable/result/FFPE_FF_MiSeq_gene.count",
              "H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/miseq/star_deseq2/result/MiSeq_FFPE_VS_FF.design", 
              "MiSeq_FF", "MiSeq_FFPE")
)

library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")
library(plyr)

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
  library(grid)
  
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

comparisonNames=names(comparisons)
comparisonName=comparisonNames[2]

tnbcgenes<-read.table("E:/sqh/Dropbox/sciences/projects/20150226_bojana_MiSeq_HiSeq/TNBC_geneid.txt", header=F, stringsAsFactors=F)$V1

protein_coding_only<-TRUE
#for(tnbcgenesonly in c(TRUE, FALSE)){
for(tnbcgenesonly in c(FALSE)){
  corrdata<-data.frame(Dataset=character(), 
                       Corr=numeric(),
                       stringsAsFactors=FALSE)
  for(comparisonName in comparisonNames){
    str(comparisonName)
    
    dataFile=comparisons[[comparisonName]][1]
    designFile=comparisons[[comparisonName]][2]
    gnames=comparisons[[comparisonName]][3:4]
    
    data<-read.table(dataFile,row.names=1, header=T, check.names=F)
    if(tnbcgenesonly){
      data<-data[data$Feature_gene_name %in% tnbcgenes,]
    }
    
    countData<-data
    
    if(protein_coding_only){
      countData<-countData[grepl("protein_coding", countData$Feature_gene_biotype),]
    }
    
    countData<-countData[,!grepl("Feature_", colnames(countData))]
    
    countData[is.na(countData)] <- 0
    countData<-round(countData)
    
    designData<-read.table(designFile, sep="\t", header=T)
    designData$Condition<-factor(designData$Condition, levels=gnames)
    
    comparisonData<-countData[,as.character(designData$Sample),drop=F]
    if(ncol(comparisonData) != nrow(designData)){
      warning(paste0("Data not matched, there are ", nrow(designData), " samples in design file ", designFile, " but ", ncol(comparisonData), " samples in data "))
      next
    }
    
    notEmptyData<-apply(comparisonData, 1, max) > 0
    comparisonData<-comparisonData[notEmptyData,]
    
    colnames(comparisonData)<-unlist(lapply(c(1:ncol(comparisonData)), function(i){paste0(designData$Paired[i], "_", colnames(comparisonData)[i])}))
    rownames(designData)<-colnames(comparisonData)
    
    ncond<-table(designData$Condition)[1]
    nlen<-nrow(designData)
    
    corr<-apply(comparisonData, 1, function(x){
      cor(x[1:ncond], x[(ncond+1):nlen], method="spearman")
    })
    
    curcorr<-data.frame(Dataset=comparisonName, Corr=corr)
    corrdata<-rbind(corrdata, curcorr)
    #   
    #   Y=data[data$Group=="Unpaired",]$OverlapRate
    #   X=data[data$Group=="Paired",]$OverlapRate
    #   pvalue<-wilcox.test(X,Y)$p.value
    #   
    #   if(pvalue < 0.001){
    #     pvaluestr<-"p<0.001"
    #   }else{
    #     pvaluestr<-sprintf("p=%.3f", pvalue)
    #   }
    #   
    #   maxRate = round(max(data$OverlapRate)*10)/10
    #   
    #   data2 <- data.frame(x = c(1, 1, 2, 2), y = c(maxRate,maxRate+0.02,maxRate+0.02,maxRate))
    #   
    #   png(file=paste0(file, ".png"), height=3000, width=4000, res=300)
    #   g<-ggplot(data, aes(x=Group, y=OverlapRate), fill=Group)+ geom_violin() + 
    #     geom_path(data=data2, aes(x=x,y=y))+
    #     annotate("text",x=1.5,y=maxRate+0.035,label=pvaluestr)+
    #     ggtitle(basename(file))
    #   print(g)
    #   dev.off()
    
  }
  
  corrdata$Dataset<-factor(corrdata$Dataset)
  
  filename=ifelse(tnbcgenesonly,"FFPE_FF.tnbc","FFPE_FF")
  filename = ifelse(protein_coding_only,paste0(filename,"_ProteinCoding"), filename )
  filename=paste0(filename, ".corr.png")
  
  if(!file.exists(filename)){
    png(filename=filename, width=3000, height=3000, res=300)
    p1<-ggplot(corrdata, aes(Corr, fill = Dataset)) + geom_density(alpha = 0.2)
    p2<-ggplot(corrdata, aes(y = Corr, x = Dataset, fill = Dataset)) + geom_violin(alpha = 0.2)
    multiplot(p1, p2, cols=1)
    dev.off()
  }
  wilcox.test(corrdata[corrdata$Dataset=="HiSeq","Corr"], corrdata[corrdata$Dataset=="MiSeq","Corr"])
  
  cuts<-unlist(lapply(corrdata$Corr, cut, c(-Inf, 0.8, 0.85, 0.9, 0.95, Inf), labels=0:4))
  corrdata$Category<-revalue(cuts, c("0"="<0.8", "1"="0.8~0.85", "2"="0.85~0.9", "3"="0.9~0.95", "4"=">0.95"))
  
  tb<-table(corrdata$Dataset, corrdata$Category)
  tb[1,]<-100*tb[1,]/sum(tb[1,])
  tb[2,]<-100*tb[2,]/sum(tb[2,])
  
  filename=ifelse(tnbcgenesonly,"FFPE_FF.tnbc.percentage.csv","FFPE_FF.percentage.csv")
  write.csv(tb, file=filename)
}