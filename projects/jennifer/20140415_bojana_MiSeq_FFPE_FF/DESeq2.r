library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")

setwd("H:/shengquanhu/projects/Jennifer/20140415_bojana_MiSeq_FFPE_FF")

countData<-read.csv("IG-20140404-count_table-gene_symbol_FF_FFPE_matched.csv",row.names=1, header=T, check.names=F)
countData<-countData[,c(5,8,17,18,1:4,6,7,9:16)]
countData<-countData[apply(countData, 1, max) > 0,]
colnames(countData)<-unlist(lapply(colnames(countData), function(x){
  substr(x, 1, nchar(x) - 6)
}))

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

colData=read.csv("MiSeq_Files_IDs_Pietenpol_lab_4_4_2014_forTiger_update.csv",row.names=1, header=T) 
colData<-colData[,c(1,2,7)]
colnames(colData)<-c("PatientID", "SampleType", "Age")

countData<-countData[,rownames(colData)]

colors<-as.matrix(data.frame(SampleType=rep(c("red", "blue"), 9), Age=c(rep(c("red"), 10),rep(c("blue"), 8))))
    
dds=DESeqDataSetFromMatrix(countData = countData,
                           colData = colData,
                           design = ~ PatientID + Age + SampleType)


dds <- DESeq(dds)
    res<-results(dds,cooksCutoff=FALSE)
    
    cat("DESeq2 finished.\n")
    
    select<- (!is.na(res$padj)) & (res$padj<0.05) & ((res$log2FoldChange >= 1) | (res$log2FoldChange <= -1))
    
    if(length(indecies) > 0){
      tbb<-cbind(data[,indecies,drop=F], pairCountData, res)
    }else{
      tbb<-cbind(pairCountData, res)
    }
    tbbselect<-tbb[select,,drop=F]
    
    tbb<-tbb[order(tbb$padj),,drop=F]
    write.csv(as.data.frame(tbb),paste0(pairname, "_DESeq2.csv"))
    
    tbbselect<-tbbselect[order(tbbselect$padj),,drop=F]
    write.csv(as.data.frame(tbbselect),paste0(pairname, "_DESeq2_sig.csv"))
    #some basic graph
    dds=DESeqDataSetFromMatrix(countData = pairCountData,
                               colData = pairColData,
                               design = ~1)
    
    colnames(dds)<-colnames(pairCountData)
    
    #draw density graph
    rldmatrix<-as.matrix(log2(counts(dds,normalized=FALSE) + 1))
    rsdata<-melt(rldmatrix)
    colnames(rsdata)<-c("Gene", "Sample", "log2Count")
    png(filename=paste0(pairname, "_DESeq2-log2-density.png"), width=4000, height=3000, res=300)
    g<-ggplot(rsdata) + geom_density(aes(x=log2Count, colour=Sample)) + xlab("DESeq2 log2 transformed count")
    print(g)
    dev.off()
    
    #varianceStabilizingTransformation
    vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
    assayvsd<-assay(vsd)
    write.csv(assayvsd, file=paste0(pairname, "_DESeq2-vsd.csv"))
    
    vsdiqr<-apply(assayvsd, 1, IQR)
    assayvsd<-assayvsd[order(vsdiqr, decreasing=T),]
    
    rldmatrix=as.matrix(assayvsd)
    
    #draw pca graph
    png(filename=paste0(pairname, "_DESeq2-vsd-pca.png"), width=3000, height=3000, res=300)
    pca<-prcomp(t(rldmatrix))
    supca<-summary(pca)$importance
    pcadata<-data.frame(pca$x)
    pcalabs=paste0(colnames(pcadata), "(", round(supca[2,] * 100), "%)")
    g <- ggplot(pcadata, aes(x=PC1, y=PC2, label=row.names(pcadata))) +
      geom_text(vjust=-0.6, size=4) +
      geom_point(col=pairColors, size=4) +
      scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) +
      scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) +
      geom_hline(aes(0), size=.2) +
      geom_vline(aes(0), size=.2) +
      xlab(pcalabs[1]) + ylab(pcalabs[2])
    print(g)
    dev.off()
    
    #draw heatmap
    rldselect<-rldmatrix[1:500,,drop=F]
    if(nrow(rldselect) > 2){
      png(filename=paste0(pairname, "_DESeq2-vsd-heatmap.png"), width=3000, height =3000, res=300)
      heatmap3(rldselect, col = hmcols, ColSideColors = pairColors, margins=c(12,5), scale="r", dist=dist, lab
               Row="",
               legendfun=function() showLegend(legend=paste0("Group ", gnames),col=c("red","blue"),cex=
                                                 1.5,x="center"))
      dev.off()
    }
    
  }
  