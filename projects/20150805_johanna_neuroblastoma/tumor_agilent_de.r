setwd("H:/shengquanhu/projects/Jennifer/20150805_johanna_neuroblastoma")

library(limma)  

data<-read.csv("GSE49710_Probe_Expression.csv", header=T, row.names=1)

candidates<-c("MYC","MYCL1", "MYCN", "TFPI2")

#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#F0E442", "red", "#0072B2", "#D55E00", "#CC79A7",  "#009E73", "blue")

drawpoints<-function(ddata, filename, title){
  ddata$Gene<-rownames(ddata)
  meltdata<-melt(ddata, id="Gene")
  colnames(meltdata)<-c("Gene","Sample","Expression")
  png(file=filename, height=1300, width=max(2000, 50 * nsample), res=300)
  g<-ggplot(meltdata, aes(x=Sample, y=Expression))+ 
    geom_point(aes(colour=Gene), size=3) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10, hjust=0, face="bold", color="black"),
          axis.text.y = element_text(size=10, face="bold", color="black"),
          plot.title=element_text(face="bold", size=18))
  print(g)
  dev.off()
}

pngfile<-"GSE49710_ByMYC.png"
if(!file.exists(pngfile)){
  probedata<-data[data$SYMBOL %in% candidates, ]
  probedata$Mean<-apply(probedata[,-1],1,mean)
  
  gene<-candidates[1]
  
  fdata<-NULL
  for(gene in candidates){
    pdata<-probedata[probedata$SYMBOL==gene,]  
    maxmean = max(pdata$Mean)
    kdata<-pdata[pdata$Mean == maxmean,]
    if(is.null(fdata)){
      fdata<-kdata
    }else{
      fdata<-rbind(fdata, kdata)
    }
  }
  
  rownames(fdata)<-fdata$SYMBOL
  fdata<-fdata[,2:(ncol(fdata)-1)]
  nsample<-ncol(fdata)
  
  #order by MYC
  ddata<-fdata[,order(fdata["MYC",])]
  drawpoints(ddata, pngfile, "GSE16476 Gene Expression Values Order by MYC")
  
  #order by MYCN
  ddata<-fdata[,order(fdata["MYCN",])]
  drawpoints(ddata, "GSE49710_ByMYCN.png", "GSE49710 Gene Expression Values Order by MYCN")
  
  #order by MYCN
  ddata<-fdata[,order(fdata["TFPI2",])]
  drawpoints(ddata, "GSE49710_ByTFPI2.png", "GSE49710 Gene Expression Values Order by TFPI2")
  
  distances<-unlist(apply(fdata, 2, function(x){
    x[1] - x[3]
  }))
  ddata<-fdata[,order(distances)]
  drawpoints(ddata, "GSE49710_ByMYC_MYCN_Distance.png", "GSE49710 Gene Expression Values Order by Distance Between MYC and MYCN")
}

canddata<-data[,-1]
nsam<-round(ncol(canddata) / 10)

defile<-"GSE49710_MYC_10percent.tsv"
if(!file.exists(defile)){
  dedata<-canddata[,order(canddata["UKv4_A_23_P215956",])]
  dedata<-dedata[,c(1:nsam, (ncol(dedata)-nsam+1):ncol(dedata))]
  
  combn<-factor(c(rep("MYC_LOW",nsam), rep("MYC_HIGH",nsam)))
  design <- model.matrix(~combn) 
  
  fit <- lmFit(dedata, design)
  efit <- eBayes(fit)
  tp<-topTable(efit, coef=2, number=nrow(dedata))
  tp<-tp[tp$adj.P.Val < 0.05,]
  tp<-tp[abs(tp$logFC) >= 1,]
  
  finaldata<-cbind(data.frame(Probe=rownames(tp), Gene=data[rownames(tp),"SYMBOL"]), tp[,c(1,4,5)], dedata[rownames(tp),])
  write.table(finaldata, file=defile, row.names=F, sep="\t")
}


defile<-"GSE49710_MYC_MYCN_DISTANCE_10percent.tsv"
if(!file.exists(defile)){
  pdata<-canddata[c("UKv4_A_23_P215956", "UKv4_A_24_P94402"),]
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
  
  finaldata<-cbind(data.frame(Probe=rownames(tp), Gene=data[rownames(tp),"SYMBOL"]), tp[,c(1,4,5)], dedata[rownames(tp),])
  write.table(finaldata, file=defile, row.names=F, sep="\t")
}