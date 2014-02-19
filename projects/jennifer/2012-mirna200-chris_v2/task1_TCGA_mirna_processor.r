source("e:/sqh/programs/r/mirna/chris_v2/common.r")

selectLowHighSamples<-function(targetRankOrdered, exprTarget, slimnames, slimfile, slimvaluefile, clusterfile, hmcols, width=3000, height=2000 ){
  tumornames<-getTumorNames(targetRankOrdered)
  samplenames<-getSampleTumors(targetRankOrdered)
  
  lowest<-c()
  highest<-c()
  percent<-0.05
  for(tumor in slimnames){
    #tumor<-"brca"
    tumorData<-targetRankOrdered[,samplenames==tumor]
    allcount<-ncol(tumorData)
    count<-round(allcount * percent)
    lowest<-c(lowest, colnames(tumorData)[1:count])
    highest<-c(highest, colnames(tumorData)[(allcount-count+1):allcount])
  }
  
  lowestdata<-targetRankOrdered[, colnames(targetRankOrdered) %in% lowest]
  highestdata<-targetRankOrdered[, colnames(targetRankOrdered) %in% highest]
  selecteddata<-cbind(lowestdata, highestdata)
  
  lowhigh<-data.frame(low=lowest, high=highest)
  write.table(lowhigh, file=slimfile, row.names=F, sep="\t")
  saveData(selecteddata, slimvaluefile)
  
  #draw heatmap by target miRNAs on high/low expressed groups.
  exprSelected<-cbind(exprTarget[,colnames(exprTarget) %in% lowest], exprTarget[,colnames(exprTarget) %in% highest])
  drawHeatmapPlusTumorAndExpressionGroup(exprSelected, clusterfile, heatmapcolor= hmcols,  width=width, height=height, cexRow=0.9)
}

#read mirna data
mirnaCount<-loadOriginalData("mirnaseq", "count")
brcaIndex<-grep("brca_", colnames(mirnaCount), value=FALSE)
brcaCount<-mirnaCount[,brcaIndex]
saveData(brcaCount, "TCGA_mirnaseq_20130521_brca_count.csv")

#keep only observed in at least 90% samples
samplecount<-ncol(brcaCount) * 0.9
keep<-rowSums(brcaCount > 0) >= samplecount

#draw heatmap by filtered miRNA
brca90Count<-brcaCount[keep,]
saveData(brca90Count, "TCGA_mirnaseq_20130521_brca_percentage90_count.csv")

brca90rank<-apply(brca90Count, 2, function(x) rank(x,ties.method="average"))
saveData(brca90rank, "TCGA_mirnaseq_20130521_brca_percentage90_rank.csv")

targetBrca90Rank<-brca90rank[row.names(brca90rank) %in% targetnames,]
saveData(targetBrca90Rank, "TCGA_mirnaseq_20130521_brca_percentage90_target_rank.csv")

targetBrca90RankOfSample<-apply(targetBrca90Rank, 1, function(x) rank(x,ties.method="average"))
saveData(targetBrca90Rank, "TCGA_mirnaseq_20130521_brca_percentage90_target_rank_sample.csv")

