source("e:/sqh/programs/r/mirna/chris/common.r")

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

targetCount<-mirnaCount[rownames(mirnaCount) %in% targetnames,]
saveData(targetCount, "TCGA_mirnaseq_20130521_count_no_ov_min100.target.csv")

#sample table
samplefile<-"TCGA_20130521_samples_no_ov_min100.csv"
if(!file.exists(samplefile)){
  Tumor<-getSampleTumors(mirnaCount)
  tumortable<-table(Tumor)
  write.csv(tumortable, samplefile, row.names=F)
}

#using DESeq to generate variance stabilisation transformed data for heatmap
tumor<-getSampleTumors(mirnaCount)
design<-data.frame(row.names=colnames(mirnaCount), tumor=tumor)
cds<-newCountDataSet(mirnaCount, design)
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds, method="blind", fitType="local")
vsd<-varianceStabilizingTransformation(cds)
expr<-exprs(vsd)

#keep only observed in at least 90% samples
samplecount<-ncol(mirnaCount) * 0.9
keep<-rowSums(mirnaCount > 0) >= samplecount

#draw heatmap by filtered miRNA
expr90<-expr[keep,]
expr90png<-"TCGA_mirnaseq_20130521_countexpr_no_ov_min100.t.90.DESeq.heatmap.png"
if(!file.exists(expr90png)){
  drawHeatmapPlusTumor(expr90, expr90png, heatmapcolor= hmcols, width=4000, height=2000, hidegenes=TRUE)
}

#draw heatmap by target miRNAs. lgg and kirc show different pattern with others
exprTarget<-expr[rownames(expr) %in% targetnames,]
exprTargetPng<-"TCGA_mirnaseq_20130521_countexpr_no_ov_min100.t.90.target.DESeq.heatmap.png"
if(!file.exists(exprTargetPng)){
  drawHeatmapPlusTumor(exprTarget, exprTargetPng, heatmapcolor= hmcols, cexRow=0.8, width=4000, height=1000)
  gc<-getTumorColors(exprTarget)
  x<-rep(1, length(gc))
  png("TCGA_mirnaseq_20130521_countexpr_no_ov_min100.t.90.target.DESeq.heatmap.color.png", width=500, height=1000, res=300)
  par(mar=c(2,3,0,1) + 0.1)
  barplot(x, col=gc, names.arg=names(gc), horiz=T, axes=F, cex.names=0.9, las=1)
  dev.off()
}

#filtered by expr value
prefix<-getFile("mirnaseq", "countexpr", "")

expr90rank<-apply(expr90, 2, function(x) rank(x,ties.method="average"))
targetExpr90Rank<-expr90rank[row.names(expr90rank) %in% targetnames,]
targetExpr90RankMean<-apply(targetExpr90Rank, 2, function(x) mean(x))
targetExpr90RankOrdered<-rbind(targetExpr90Rank, targetExpr90RankMean)
targetExpr90RankOrdered<-targetExpr90RankOrdered[,order(targetExpr90RankMean)]
rownames(targetExpr90RankOrdered)[6] <- "meanOfRank"
saveData(targetExpr90RankOrdered, paste0(prefix,".t.90.rank.target.csv"))

exprTargetRanked<-exprTarget[,order(targetExpr90RankMean)]
drawHeatmapPlusTumor(exprTargetRanked, "TCGA_mirnaseq_20130521_countexpr_no_ov_min100.t.90.target.ranked.DESeq.heatmap.png", heatmapcolor= hmcols, width=4000, height=1000, cexRow=0.8,Colv=NA)

#select low/high expressed groups
tumornames<-getTumorNames(targetExpr90RankOrdered)
slimnames<-tumornames
selectLowHighSamples(targetExpr90RankOrdered, 
                     exprTarget, 
                     slimnames,
                     paste0(prefix,".t.90.rank.target.lowhigh.sample.tsv"),
                     paste0(prefix,".t.90.rank.target.lowhigh.rank.csv"),
                     paste0(prefix,".t.90.rank.target.lowhigh.DESeq.heatmap.png"),
                     hmcols)

#remove lgg, kirc based on cluster result of target miRNAs.
slimnames<-tumornames[-which(tumornames=="kirc" | tumornames=="lgg")]
selectLowHighSamples(targetExpr90RankOrdered, 
                     exprTarget, 
                     slimnames,
                     paste0(prefix,".t.90.rank.target.lowhigh_no_igg_kirc.sample.tsv"),
                     paste0(prefix,".t.90.rank.target.lowhigh_no_igg_kirc.rank.csv"),
                     paste0(prefix,".t.90.rank.target.lowhigh_no_igg_kirc.DESeq.heatmap.png"),
                     hmcols,
                     height=2000)

######################using count and rpkm for rank get exact same result#############################
# vtypes<-c("count", "rpkm")
# for(v in vtypes){
#   vtype<-v
#   #vtype<-"count"
#   
#   prefix<-getFile("mirnaseq", vtype, "")
#   
#   #load original data
#   mirna<-loadOriginalData("mirnaseq", vtype)
#   
#   #filter target miRNAs
#   mirnaTargetFile<-paste0(prefix,".t.target.csv")
#   mirnaTarget<-mirna[row.names(mirna) %in% targetnames,]
#   saveData(mirnaTarget,mirnaTargetFile)
#   
#   #filter miRNAs based on zero count, we asked for at least the mirna detected by at least 90% samples
#   mirna90File<-paste0(prefix,".t.90.csv")
#   mirna90<-mirna[keep,]
#   saveData(mirna90,mirna90File)
#   
#   #get rank of target miRNAs
#   targetFile<-paste0(prefix,".t.90.rank.target.csv")
#   mirna90rank<-apply(mirna90, 2, function(x) rank(x,ties.method="average"))
#   targetRank<-mirna90rank[row.names(mirna90rank) %in% targetnames,]
#   targetRankMean<-apply(targetRank, 2, function(x) mean(x))
#   targetRankOrdered<-rbind(targetRank, targetRankMean)
#   targetRankOrdered<-targetRankOrdered[,order(targetRankMean)]
#   rownames(targetRankOrdered)[6] <- "meanOfRank"
#   saveData(targetRankOrdered, targetFile)
#   
#   targetValueFile<-paste0(prefix,".t.90.rank.target.value.csv")
#   targetValue<-mirna90[row.names(mirna90) %in% targetnames,order(targetRankMean)]
#   saveData(targetValue, targetValueFile)
#   
#   #select low/high expressed groups
#   tumornames<-getTumorNames(targetRankOrdered)
#   slimnames<-tumornames
#   selectLowHighSamples(targetRankOrdered, 
#                        exprTarget, 
#                        slimnames,
#                        paste0(prefix,".t.90.rank.target.lowhigh.tsv"),
#                        paste0(prefix,".t.90.rank.target.lowhigh.rank.csv"),
#                        paste0(prefix,".t.90.rank.target.lowhigh.DESeq.heatmap.png"),
#                        hmcols)
#   
#   #remove lgg, kirc based on cluster result of target miRNAs.
#   slimnames<-tumornames[-which(tumornames=="kirc" | tumornames=="lgg")]
#   selectLowHighSamples(targetRankOrdered, 
#                        exprTarget, 
#                        slimnames,
#                        paste0(prefix,".t.90.rank.target.lowhigh_no_igg_kirc.tsv"),
#                        paste0(prefix,".t.90.rank.target.lowhigh_no_igg_kirc.rank.csv"),
#                        paste0(prefix,".t.90.rank.target.lowhigh_no_igg_kirc.DESeq.heatmap.png"),
#                        hmcols)
# }
