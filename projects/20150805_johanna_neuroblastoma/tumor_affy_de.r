setwd("H:/shengquanhu/projects/Jennifer/20150805_johanna_neuroblastoma")

library(DupChecker)
#Microarray
#geoDownload(c("GSE16476"))

#RNASeq
#geoDownload(c("GSE49711"))

#Microarray
#geoDownload(c("GSE49710"))

library(affy)   # Affymetrix pre-processing
library(limma)  # two-color pre-processing; differential
library(hgu133plus2.db)
library(heatmap3)

evaluesfile<-"GSE16476.csv"
if(file.exists(evaluesfile)){
  evalues<-read.csv(evaluesfile, row.names=1)
}else{
  ## RMA normalization
  celfiles <- list.files("GSE16476",".gz$")
  eset <- justRMA(filenames=celfiles, celfile.path="GSE16476")
  
  evalues<-exprs(eset)
  colnames(evalues)<-substr(colnames(evalues), 1, 9)
  write.csv(evalues, file=evaluesfile)
}

es<-select(hgu133plus2.db, rownames(evalues), "SYMBOL", "PROBEID")
es<-es[!is.na(es$SYMBOL),]

uniqueProbe<-unique(es$PROBEID)
uniqueProtein<-unlist(lapply(uniqueProbe, function(x){
  ues<-unlist(es[es$PROBEID == x,]$SYMBOL)
  paste(ues,collapse=";")
}))

es<-data.frame(PROBEID=uniqueProbe, SYMBOL=uniqueProtein)
rownames(es)<-es$PROBEID

candidates<-c("MYC","MYCL", "MYCN", "TFPI2")

#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#F0E442", "red", "#0072B2", "#D55E00", "#CC79A7",  "#009E73", "blue")

drawpoints<-function(ddata, filename, title){
  ddata$Gene<-rownames(ddata)
  meltdata<-melt(ddata, id="Gene")
  colnames(meltdata)<-c("Gene","Sample","RMA")
  png(file=filename, height=1300, width=max(2000, 50 * nsample), res=300)
  g<-ggplot(meltdata, aes(x=Sample, y=RMA))+ 
    geom_point(aes(colour=Gene), size=3) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10, hjust=0, face="bold", color="black"),
          axis.text.y = element_text(size=10, face="bold", color="black"),
          plot.title=element_text(face="bold", size=18))
  print(g)
  dev.off()
}

pngfile<-"GSE16476_ByMYC.png"
if(!file.exists(pngfile)){
  probes<-es[es$SYMBOL %in% candidates, "PROBEID"]
  
  probedata<-evalues[rownames(evalues) %in% probes,]
  probedata$Mean<-apply(probedata,1,mean)
  probedata$Gene<-as.character(es[rownames(probedata),"SYMBOL"])
  
  gene<-candidates[1]
  
  fdata<-NULL
  for(gene in candidates){
    pdata<-probedata[probedata$Gene==gene,]  
    maxmean = max(pdata$Mean)
    kdata<-pdata[pdata$Mean == maxmean,]
    if(is.null(fdata)){
      fdata<-kdata
    }else{
      fdata<-rbind(fdata, kdata)
    }
  }
  
  rownames(fdata)<-fdata$Gene
  fdata<-fdata[,1:(ncol(fdata)-2)]
  nsample<-ncol(fdata)
  
  #order by MYC
  ddata<-fdata[,order(fdata["MYC",])]
  drawpoints(ddata, pngfile, "GSE16476 Gene Expression Values Order by MYC")
  
  #order by MYCN
  ddata<-fdata[,order(fdata["MYCN",])]
  drawpoints(ddata, "GSE16476_ByMYCN.png", "GSE16476 Gene Expression Values Order by MYCN")
  
  #order by MYCN
  ddata<-fdata[,order(fdata["TFPI2",])]
  drawpoints(ddata, "GSE16476_ByTFPI2.png", "GSE16476 Gene Expression Values Order by TFPI2")
  
  distances<-unlist(apply(fdata, 2, function(x){
    x[1] - x[3]
  }))
  ddata<-fdata[,order(distances)]
  drawpoints(ddata, "GSE16476_ByMYC_MYCN_Distance.png", "GSE16476 Gene Expression Values Order by Distance Between MYC and MYCN")
}

canddata<-evalues[rownames(evalues) %in% es$PROBEID,]
nsam<-round(ncol(canddata) / 10)

defile<-"GSE16476_MYC_10percent.tsv"
if(!file.exists(defile)){
  dedata<-canddata[,order(canddata["202431_s_at",])]
  dedata<-dedata[,c(1:nsam, (ncol(dedata)-nsam+1):ncol(dedata))]
  
  combn<-factor(c(rep("MYC_LOW",nsam), rep("MYC_HIGH",nsam)))
  design <- model.matrix(~combn) 
  
  fit <- lmFit(dedata, design)
  efit <- eBayes(fit)
  tp<-topTable(efit, coef=2, number=nrow(dedata))
  tp<-tp[tp$adj.P.Val < 0.05,]
  tp<-tp[abs(tp$logFC) >= 1,]
  
  finaldata<-cbind(data.frame(Probe=rownames(tp), Gene=es[rownames(tp),"SYMBOL"]), tp[,c(1,4,5)], dedata[rownames(tp),])
  write.table(finaldata, file=defile, row.names=F, sep="\t")
}


defile<-"GSE16476_MYC_MYCN_DISTANCE_10percent.tsv"
if(!file.exists(defile)){
  pdata<-canddata[c("202431_s_at", "209757_s_at"),]
  distances<-unlist(apply(pdata, 2, function(x){
    x[1] - x[2]
  }))
  dedata<-canddata[,order(distances)]
  dedata<-dedata[,c(1:nsam, (ncol(dedata)-nsam+1):ncol(dedata))]
  
  combn<-factor(c(rep("MYC_LOW",nsam), rep("MYC_HIGH",nsam)))
  design <- model.matrix(~combn) 
  
  fit <- lmFit(dedata, design)
  efit <- eBayes(fit)
  tp<-topTable(efit, coef=2, number=nrow(dedata))
  tp<-tp[tp$adj.P.Val < 0.05,]
  tp<-tp[abs(tp$logFC) >= 1,]
  
  finaldata<-cbind(data.frame(Probe=rownames(tp), Gene=es[rownames(tp),"SYMBOL"]), tp[,c(1,4,5)], dedata[rownames(tp),])
  write.table(finaldata, file=defile, row.names=F, sep="\t")
}