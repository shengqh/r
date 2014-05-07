library("sfsmisc")
library("WGCNA")
options(stringsAsFactors = FALSE);

setwd("H:/shengquanhu/projects/chenxi/20130724-chenxi_colon_cancer");

#load the candidate probe~gene list
genes<-read.csv("gene_selection_stage23.csv", header=T, row.names=1)

dataname<-"training.rma.combat"

#load expression value
load(paste0(dataname, ".RData"))

#keep candidate genes only
wdata<-t(batchExpr[row.names(genes),])
colnames(wdata)<-genes[colnames(wdata),"anno.geneSymbol"]
wname<-paste0(dataname, ".survialgenes")
save(wdata, file=paste0(wname, ".RData"))

rm(dataname)
rm(batchExpr)

#find power
powers = c(1:10)
sft = pickSoftThreshold(wdata, powerVector = powers, verbose = 5)[[2]]

png(paste0(wname, ".findpower.png"), width=4000, height=3000,res=300)

par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft[,1], -sign(sft[,3])*sft[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft[,1], -sign(sft[,3])*sft[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft[,1], sft[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft[,1], sft[,5], labels=powers, cex=cex1,col="red")

dev.off()

#based on the image, 7 is the power
powertolerance=7

tomfile<-paste0(wname, ".TOM")
net = blockwiseModules(wdata, power = powertolerance,
                       saveTOMs = FALSE,
                       verbose = 7)

# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

colors<-unique(as.character(moduleColors))
save(colors, file=paste0(wname, ".colors.RData"))


# Plot the dendrogram and the module colors underneath
png(paste0(wname, ".cluster.png"), width=4000, height=3000,res=300)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    groupLabels = "Modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs = net$MEs;
geneTree = net$dendrograms[[1]];

save(net, file = paste0(wname, ".netauto.RData"))

#Calculate topological overlap
TOM = TOMsimilarityFromExpr(wdata, power = powertolerance);
dimnames(TOM) = list(colnames(wdata), colnames(wdata))

#get distance
dissTOM = 1-TOM

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^powertolerance

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

# Call the plot function
png(paste0(wname, ".heatmap.png"), width=4000, height=3000,res=300)
TOMplot(plotTOM, geneTree, moduleColors, main = paste0("Network heatmap plot, ", ncol(wdata), " genes"))
dev.off()

# Export to cytoscape
modColors<-c("all", unique(as.character(moduleColors)))

thresholds<-c(0.02, 0.05)
thre<-0.2


thres<-c(1:30) / 100

counts<-list()

modColor<-colors[1]
for(modColor in colors){
  select<-moduleColors==modColor
  modTOM<-TOM[select, select]
  # Export the network into edge and node list files Cytoscape can read
  modGenes<-rownames(modTOM)
  counts[[modColor]]<-unlist(lapply(thres, function(x){
    cyt = exportNetworkToCytoscape(modTOM,
                                   weighted = TRUE,
                                   threshold = x,
                                   nodeAttr = as.character(moduleColors[select]))
    nrow(cyt[["nodeData"]])
  }))
}
modulecounts<-data.frame(counts)

png(paste0(wname, ".thresholds.png"), width=4000, height=3000,res=300)
matplot(thres, modulecounts, type="l", lty=1, col=colors, xlab="Threshold", ylab="Module Count",pch=c(1:length(colors)))
matpoints(thres, modulecounts, lty=1, col=colors,pch=c(1:length(colors)))
legend("topright", legend=colors, col=colors, pch=c(1:length(colors)))
dev.off()

for(modColor in modColors){
  if(modColor=="all"){
    select<-rep(TRUE, ncol(TOM))
  }else{
    select<-moduleColors==modColor
  }
  # Select modules
  modTOM<-TOM[select, select]
  # Export the network into edge and node list files Cytoscape can read
  modGenes<-rownames(modTOM)
  
  for(thre in thresholds){
    edgeFile = paste0(wname, "-", paste(modColor, collapse="-"), "-", thre, "-edges.tsv")
    nodeFile = paste0(wname, "-", paste(modColor, collapse="-"), "-", thre, "-nodes.tsv")
    
    cyt = exportNetworkToCytoscape(modTOM,
                                   edgeFile = edgeFile,
                                   nodeFile = nodeFile,
                                   threshold = thre,
                                   nodeAttr = as.character(moduleColors[select]))
    
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