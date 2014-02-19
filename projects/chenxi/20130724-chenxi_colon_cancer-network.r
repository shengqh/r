require("WGCNA")

setwd("H:/shengquanhu/projects/chenxi/20130724-chenxi_colon_cancer");

name<-"coloncancer_training_stage23.csv.qnorm.bnorm.pregenes"
load(paste0(name, ".RData"))
load(paste0(name, ".netauto.RData"))

powertolerance<-6

datExpr=t(pregdata)

# Export to cytoscape
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = powertolerance);
dimnames(TOM) = list(rownames(pregdata), rownames(pregdata))

moduleColors[moduleColors=="grey"]<-"red"
moduleColors[moduleColors=="turquoise"]<-"green"
modColors<-unique(moduleColors)

thresholds<-c(0.02, 0.1, 0.2)
thre<-0.2
for(modColor in modColors){
  select<-moduleColors==modColor
  # Select modules
  modTOM<-TOM[select, select]
  # Export the network into edge and node list files Cytoscape can read
  modGenes<-rownames(modTOM)
  
  for(thre in thresholds){
    edgeFile = paste0(name, "-", paste(modColor, collapse="-"), "-", thre, "-edges.tsv")
    nodeFile = paste0(name, "-", paste(modColor, collapse="-"), "-", thre, "-nodes.tsv")
    
    cyt = exportNetworkToCytoscape(modTOM,
                                   edgeFile = edgeFile,
                                   nodeFile = nodeFile,
                                   weighted = TRUE,
                                   threshold = thre,
                                   nodeNames = modGenes,
                                   nodeAttr = rep(modColor, nrow(modTOM)))
    
    edge<-read.table(edgeFile, sep=" ", header=F, stringsAsFactors=F)
    write.table(edge, file=edgeFile, row.names=F, col.names=F, sep="\t", quote=F)

    node<-data.frame(read.table(nodeFile, sep=" ", header=F, stringsAsFactors=F))
    nodeName<-node$V1[-1]
    if(nrow(edge) > 1 ){
      edgeTable<-data.frame(table(edge$V1[-1]))
      rownames(edgeTable)<-edgeTable$Var1
      nodeConn<-edgeTable[nodeName,]$Freq
      nodeConn[is.na(nodeConn)]<-0
      node$V4<-c("nodeEdgeNumber", nodeConn  )
    }else{
      node$V4<-c("nodeEdgeNumber", rep(0, nrow(node)-1)  )
    }
    node<-node[c(1, order(as.numeric(node$V4[-1]), decreasing=T)+1),]
    node<-node[,c(1,3:ncol(node))]
    write.table(node, file=nodeFile, row.names=F, col.names=F, sep="\t", quote=F)
  }
}