source("e:/sqh/programs/r/mirna/chris/common.r")

xfile<-"TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.csv"
rnaseqCount<-loadRNAseqData("count", xfile)

#definition
ncoltop<-ncol(rnaseqCount) / 2
group<-c(rep("low", ncoltop), rep("high", ncoltop))
tumor<-getSampleTumors(rnaseqCount)

#using DESeq to generate variance stabilisation transformed data for heatmap
vsdfile<-"TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.vsd.csv"
if(file.exists(vsdfile)){
  expr<-loadData(vsdfile)
}else{
  design<-data.frame(row.names=colnames(rnaseqCount), tumor=tumor)
  cds<-newCountDataSet(rnaseqCount, design)
  cds<-estimateSizeFactors(cds)
  cds<-estimateDispersions(cds, method="blind", fitType="local")
  vsd<-varianceStabilizingTransformation(cds)
  expr<-exprs(vsd)
  saveData(expr,vsdfile)
}

cdsfile<-"TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.cds.RData"
if(file.exists(cdsfile)){
  load(cdsfile)
}else{
  design<-data.frame(row.names=colnames(rnaseqCount), tumor=tumor, mirna=group)
  cds<-newCountDataSet(rnaseqCount, design)
  cds<-estimateSizeFactors(cds)
  cds<-estimateDispersions(cds)
  save(cds,file=cdsfile)
}

#filter genes by occurance
keep<-rowSums(rnaseqCount > 0) >= ncoltop
cdsFilt<-cds[keep,]
print(paste0("There are ", nrow(cdsFilt), " genes observed in at least half of samples"))
  
#DE by edge
diffFile<-"TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.DESeq_glm.csv"
if(file.exists(diffFile)){
  diff<-loadData(diffFile)
}else{
  xx<-rowSums(counts(cdsFilt) == 0) / nrow(cdsFilt)

  fit1file<-paste0(diffFile, ".fit1.Rdata")
  if(file.exists(fit1file)){
    load(fit1file)
  }else{
    fitFilt1 = fitNbinomGLMs(cdsFilt, count ~ tumor + mirna)
    save(fitFilt1, file=fit1file)
  }

  fit0file<-paste0(diffFile, ".fit0.Rdata")
  if(file.exists(fit0file)){
    load(fit0file)
  }else{
    fitFilt0 = fitNbinomGLMs(cdsFilt, count ~ tumor)
    save(fitFilt0, file=fit0file)
  }
  
  PValue = nbinomGLMTest(fitFilt1, fitFilt0)
  FDR = p.adjust(pvalsFilt, method="bonferroni")
  
  diff<-cbind(fitFilt1, PValue, FDR)
  diff<-diff[!is.na(FDR),]
  diff<-diff[diff$converged,]
  diff<-diff[order(diff$FDR),]
  
  saveData(diff, diffFile)
}

folds<-c(1,1.5,2)
fdr<-0.01
for(useabs in c(TRUE, FALSE)){
  for(fold in folds){
    if(useabs){
      sigfile<-paste0("TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.DESeq_glm.sig_abs(fold)_", fold,".csv")
      sig<-diff[(abs(diff$mirnalow) >= fold) & (diff$FDR <= fdr),]
      foldchange<-"|log2FoldChange| >= "
    }else{
      sigfile<-paste0("TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.DESeq_glm.sig_fold_", fold,".csv")
      sig<-diff[(diff$mirnalow >= fold) & (diff$FDR <= fdr),]
      foldchange<-"log2FoldChange >= "
    }
    write.csv(sig, sigfile)
    exprsig<-expr[rownames(expr) %in% rownames(sig),]
    drawHeatmapPlusTumorAndExpressionGroup(exprsig, paste0(sigfile,".heatmap.png"), heatmapcolor= hmcols,  width=3000, height=3000, cexRow=0.8, hidegenes=TRUE, main=paste0(foldchange, fold, ", FDR <= ", fdr, ", genes = ", nrow(sig)))
  }
}

#######select top 200 as significant genes#################
# top<-200
# fold<-1.5
# sigfile<-paste0("TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.DESeq_glm.sig_abs(fold)_1.5_top", top, ".csv")
# sig<-diff[(abs(diff$mirnalow) >= fold) & (diff$FDR <= fdr),]
# sig<-sig[1:top,]
# foldchange<-"|log2FoldChange| >= "
# write.csv(sig, sigfile)
# exprsig<-expr[rownames(expr) %in% rownames(sig),]
# drawHeatmapPlusTumorAndExpressionGroup(exprsig, paste0(sigfile,".heatmap.png"), heatmapcolor= hmcols,  width=3000, height=3000, cexRow=0.8, hidegenes=TRUE, main=paste0(foldchange, fold, ", FDR <= ", fdr, ", genes = ", nrow(sig)))
