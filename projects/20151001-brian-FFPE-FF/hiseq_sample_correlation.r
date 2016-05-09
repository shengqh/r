
setwd("E:/sqh/Dropbox/career/papers/draft/201510-brian-FFPE-FF")  

comparisons=list(
  "HiSeq" = c("FFPE_FF_HiSeq_gene.count",
              "HiSeq_FFPE_VS_FF.design", 
              "HiSeq_FF", "HiSeq_FFPE")
)

library("ggplot2")

comparisonNames=names(comparisons)
comparisonName=comparisonNames[1]

finaldata<-data.frame(GeneType=character(),
                     Dataset=character(), 
                     Corr=numeric(),
                     stringsAsFactors=FALSE)
for(comparisonName in comparisonNames){
  str(comparisonName)
  
  dataFile=comparisons[[comparisonName]][1]
  designFile=comparisons[[comparisonName]][2]
  gnames=comparisons[[comparisonName]][3:4]
  
  data<-read.table(dataFile,row.names=1, header=T, check.names=F)
  designData<-read.table(designFile, sep="\t", header=T)
  designData$Condition<-factor(designData$Condition, levels=gnames)
  
  index<-1
  for(index in c(1:3)){
    protein_coding_only<-pcs[index]
    tnbcgenesonly<-tnbcs[index]
    name<-names[index]
    
    countData<-data

    if(tnbcgenesonly){
      countData<-countData[countData$Feature_gene_name %in% tnbcgenes,]
    }
    
    if(protein_coding_only){
      countData<-countData[grepl("protein_coding", countData$Feature_gene_biotype),]
    }
    
    countData<-countData[,!grepl("Feature_", colnames(countData))]
    
    countData[is.na(countData)] <- 0
    countData<-round(countData)
    
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
    
    curcorr<-data.frame(GeneType=name, Dataset=comparisonName, Corr=corr)
    finaldata<-rbind(finaldata, curcorr)
  }
}

corr<-finaldata[finaldata$GeneType==names[1],]
all_pvalue<-wilcox.test(corr[corr$Dataset=="HiSeq","Corr"], corr[corr$Dataset=="MiSeq","Corr"])$p.value

corr<-finaldata[finaldata$GeneType==names[2],]
coding_pvalue<-wilcox.test(corr[corr$Dataset=="HiSeq","Corr"], corr[corr$Dataset=="MiSeq","Corr"])$p.value

corr<-finaldata[finaldata$GeneType==names[3],]
tnbc_pvalue<-wilcox.test(corr[corr$Dataset=="HiSeq","Corr"], corr[corr$Dataset=="MiSeq","Corr"])$p.value

png(file="figure2_gene_correlation.png", width=3000, height=1600, res=300)
g<-ggplot(finaldata, aes(Corr, fill = Dataset)) + 
  geom_density(alpha = 0.2) + 
  facet_grid( . ~ GeneType) + 
  xlab("Spearman correlation of gene expression level between FF and FFPE samples") + 
  ylab("Density")+ 
  annotate("text", x = -0.5, y = 1.5, label = "p<0.001")
print(g)
dev.off()