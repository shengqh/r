library("heatmap3")

setwd("H:/shengquanhu/projects/chenxi/20130724-chenxi_colon_cancer");

genes<-read.csv("gene_selection_stage23.csv", header=T, row.names=1)

topgenes<-read.csv("topgenes.csv", header=T, row.names=1)

load("training.rma.combat.survialgenes.colors.RData")

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

datasets<-c("GSE14333", "GSE17536", "GSE17537", "GSE33113")
tols<-c("0.02", "0.05")

ds<-"GSE33113"
tol<-"0.05"
for(tol in tols){
  for(ds in datasets){
    dsfile<-paste0(ds, ".rma.RData")
    load(dsfile)
    wdata<-rmaExpr[row.names(genes),]
    rownames(wdata)<-genes[rownames(wdata),"anno.geneSymbol"]

    stageColors<-unlist(lapply(dataInfo$FinalStage, function(x){
      if(x == 2){
        "red"
      }else{
        "blue"
      }
    }))
    
    modColor<-"green"
    for(modColor in colors){
      nodes<-read.delim(paste0("training.rma.combat.survialgenes-", modColor, "-", tol, "-nodes.tsv"),row.names=1)
      nodeData<-wdata[rownames(nodes),]
      
      maxcolrow<-max(nrow(nodeData), ncol(nodeData))
      png(filename=paste0("training.rma.combat.survialgenes-", modColor, "-", tol, "-", ds, "-nodes.heatmap.png"), width=50 * maxcolrow, height =50 * maxcolrow, res=300)
      heatmap3(nodeData, col = hmcols, ColSideColors = stageColors, ColSideLabs="Stage", dist=dist,  margins=c(5,5), scale="r", legendfun=function() showLegend(legend=c("Stage2", "Stage3"),col=c("red","blue"),cex=1.5,x="center"))
      dev.off()
    }
    
    topdata<-wdata[row.names(wdata) %in% topgenes$gene,]
    maxcolrow<-max(nrow(topdata), ncol(topdata))
    png(filename=paste0("training.rma.combat.survialgenes-topgenes-", tol,"-", ds,  "-nodes.heatmap.png"), width=50 * maxcolrow, height =50 * maxcolrow, res=300)
    
    #heatmap(topdata)
    heatmap3(topdata, col = hmcols, ColSideColors = stageColors, ColSideLabs="Stage", dist=dist, margins=c(5,5), scale="r", legendfun=function() showLegend(legend=c("Stage2", "Stage3"),col=c("red","blue"),cex=1.5,x="center"))
    dev.off()
  }  
}
