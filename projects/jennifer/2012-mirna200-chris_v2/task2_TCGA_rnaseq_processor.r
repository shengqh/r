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

#using edgeR to filter genes, keep genes with at least two counts per million(CPM) in at least half of the samples
cpmx<-cpm(rnaseqCount)
keep<-rowSums(cpmx > 2) >= ncoltop
countdata<-rnaseqCount[keep,]
  
#DE by edge
diffFile<-"TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.edgeR_glm.csv"
if(file.exists(diffFile)){
  diff<-loadData(diffFile)
}else{
  xx<-rowSums(countdata == 0) / nrow(countdata)

  lrtfile<-paste0(diffFile, ".lrt.Rdata")
  if(file.exists(lrtfile)){
    load(lrtfile)
  }else{
    fitfile<-paste0(diffFile, ".fit.Rdata")
    if(file.exists(fitfile)){
      load(fitfile)
    }else{
      yfile<-paste0(diffFile, ".y.Rdata")
      if(file.exists(yfile)){
        load(yfile)
      }else{
        y <- DGEList(counts=countdata,group=group)
        design <- model.matrix(~group+tumor)
        y <- estimateGLMCommonDisp(y,design)
        y <- estimateGLMTrendedDisp(y,design)
        y <- estimateGLMTagwiseDisp(y,design)
        save(y, file=yfile)
      }
      fit <- glmFit(y,design)
      save(fit, file=fitfile)
    }
    lrt <- glmLRT(fit,coef=2)
    save(lrt, file=lrtfile)
  }

  diff<-topTags(lrt, n=nrow(countdata))$table
  diff$zeropercentage<-xx;
  diff<-diff[order(diff$FDR), ]
  saveData(diff, diffFile)
}

folds<-c(1,1.5,2)
for(useabs in c(TRUE, FALSE)){
  for(fold in folds){
    if(useabs){
      sigfile<-paste0("TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.edgeR_glm.sig_abs(fold)_", fold,".csv")
      sig<-diff[(abs(diff$logFC) >= fold) & (diff$FDR <= 0.01),]
      foldchange<-"|log2FoldChange| >= "
    }else{
      sigfile<-paste0("TCGA_rnaseqv2_20130521_count_no_ov_min100_no_igg_kirc_miRNAlowhigh.edgeR_glm.sig_fold_", fold,".csv")
      sig<-diff[(diff$logFC >= fold) & (diff$FDR <= 0.01),]
      foldchange<-"log2FoldChange >= "
    }
    write.csv(sig, sigfile)
    
    exprsig<-expr[rownames(expr) %in% rownames(sig),]
    drawHeatmapPlusTumorAndExpressionGroup(exprsig, paste0(sigfile,".heatmap.png"), heatmapcolor= hmcols,  width=3000, height=3000, cexRow=0.8, hidegenes=TRUE, main=paste0(foldchange, fold, ", FDR <= 0.01, genes = ", nrow(sig)))
  }
}
