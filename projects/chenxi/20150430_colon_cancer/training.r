library("affy");
library("preprocessCore");
library("sfsmisc")
library("genefilter")
library("sva")
library("hgu133plus2.db")
library("stringr")
library("heatmap.plus")

#function used to draw headmapplus
drawHeatmapPlus<-function(x, filename, samplecolor, heatmapcolor, width=4000, height=3000, cexRow=1.0, cexCol=1.0, main = NULL, xlab = NULL, ylab = NULL, Colv=NULL){
  if(filename != ""){
    png(filename, width=width, height=height,res=300)
  }
  par(mar=c(12, 2, 2, 2))
  heatmap.plus(x, ColSideColors = samplecolor, col = heatmapcolor,cexRow=cexRow, cexCol=cexCol, main=main, xlab=xlab, ylab=ylab, Colv=Colv)
  if(filename != ""){
    dev.off()
  }
}

#build color table of batch
getBatchColors<-function(batch){
  batchnames<-unique(batch)
  colors<-rainbow(length(batchnames))
  names(colors)<-batchnames  
  return (colors)
}

#get sample color based on batch
getSampleColors<-function(batch){
  colortable<-getBatchColors(batch)
  return (colortable[batch])
}

#draw by tumor only
drawHeatmapPlusBatch<-function(x, batch, filename="", heatmapcolor=NA, width=4000, height=3000, cexRow=1.0, cexCol=1.0, hidesamples=TRUE, hidegenes=TRUE, main = NULL, xlab = NULL, ylab = NULL, Colv = NULL){
  if(is.na(heatmapcolor)){
    heatmapcolor <- colorRampPalette(c("green", "black", "red"))(256)
  }
  
  d<-as.matrix(x)
  groupcols<-getSampleColors(batch)
  clab  <- matrix(c(rep("white", ncol(d)), groupcols), ncol=2, byrow=FALSE)
  colnames(clab)<-c("", "Batch")
  if(hidesamples){
    colnames(d)<-c(rep("", ncol(d)))
  }
  if(hidegenes){
    rownames(d)<-c(rep("", nrow(d)))
  }
  drawHeatmapPlus(d, filename, samplecolor=clab, heatmapcolor= heatmapcolor, width, height, cexRow, cexCol,main=main, xlab=xlab, ylab=ylab,Colv=Colv)
}

parseGenes<-function(rfile){
  rfileR<-paste0(rfile,".RData")
  
  gfile<-paste0(rfile,".genes")
  gfileR<-paste0(gfile, ".RData")
  
  pregfile<-paste0(rfile,".pregenes")
  pregfileR<-paste0(pregfile, ".RData")
  if(!file.exists(gfileR)){
    bdata<-local(get(load(rfileR)))
    
    cat("filter gene probe ...\n")
    
    nm<-rownames(bdata)
    isGeneProbe<-!(nm %in% grep("^A", nm, value=TRUE))
    pdata<-data.frame(subset(bdata,isGeneProbe))
    colcount<-ncol(pdata)
    #rm(bdata)
    
    pdata$hgu133plus2<-unlist(mget(rownames(pdata),envir=hgu133plus2SYMBOL))
    pdata$hgu133plus2map<-genemap[rownames(pdata),1]
    
    cat("find the best probe representing the each gene ...\n")
    
    predefined<-read.csv("final_mapping.csv", header=T, row.names=2)
    
    cat("Predefined probe = ", nrow(predefined), "\n")
    predefined$hgu133plus2<-unlist(mget(as.character(predefined$Affy.ID),envir=hgu133plus2SYMBOL))
    predefined$hgu133plus2map<-genemap[rownames(predefined),1]
    
    matched<-!is.na(predefined$hgu133plus2) & (predefined$hgu133plus2 == predefined$Nano.ID)
    unmatched<-predefined[!matched,]
    write.csv(unmatched, file="Unmatched.csv", row.names=F)
    
    #get predefined probe and gene
    predefineddata<-pdata[rownames(pdata) %in% rownames(predefined), 1:colcount]
    rownames(predefineddata)<-predefined[rownames(predefineddata), "Nano.ID"]
    
    #get postdefined probe and gene
    removegenes<-apply(predefined, 1, function(x){
      ret<-c(x["Nano.ID"], x["hgu133plus2"])
      mapr<-str_trim(unlist(strsplit(x["hgu133plus2map"], "///")))
      return (c(ret, mapr))
    })
    rgenes<-sort(unique(unlist(removegenes)))
    
    #x<-pdata["91682_at",,drop=F]
    isremove<-apply(pdata, 1, function(x){
      ret<-c(x["hgu133plus2"])
      mapr<-str_trim(unlist(strsplit(as.character(x["hgu133plus2map"]), "///")))
      curgenes<-unique(unlist(c(ret, mapr)))
      curgenes<-curgenes[!is.na(curgenes) & (curgenes != "---")]
      if(any(curgenes %in% rgenes)){
        return (TRUE)
      }else{
        return (FALSE)
      }
    })
    postdefineddata<-pdata[!isremove,]
    postdefineddata<-postdefineddata[!is.na(postdefineddata$hgu133plus2),]
    
    probes<-rownames(postdefineddata)
    
    arrayIQR<-apply(postdefineddata[,1:colcount],1,IQR)
    
    largestProbeName<-findLargest(probes, arrayIQR, "hgu133plus2")
    
    cat("present each gene by largest IQR probe ...\n")
    
    isLargestProbe<-probes %in% largestProbeName
    geneData<-subset(postdefineddata,isLargestProbe)
    #rm(postdefineddata)
    rownames(geneData)<-geneData$hgu133plus2
    
    gdata<-rbind(predefineddata, geneData[,1:colcount])
    save(gdata,file=gfileR)
    
    pregdata<-rbind(predefineddata, geneData[rownames(gdata)=="EPGN",1:colcount], geneData[rownames(gdata)=="ND4",1:colcount])
    save(pregdata, file=pregfileR)
  }
  
  pregheatfile<-paste0(pregfile, ".heatmap.png")
  if(!file.exists(pregheatfile)){
    load(pregfileR)
    
    batchdefinitionfile<-paste0(datafile,".batchdefinition");
    batchdata<-read.csv(batchdefinitionfile, header=T)
    batch<-batchdata$Batch
    
    drawHeatmapPlusBatch(pregdata, batch, file=pregheatfile)
    #drawHeatmapPlusBatch(pregdata[,1:100], batch[1:100])
  }
  
  return(c(gfile, gfileR, pregfile, pregfileR))
}

drawBoxplot<-function(data, file="", width=4000, height=2000){
  if(file != ""){
    png(file, width=width, height=height,res=300)
  }
  par(mar=c(12, 2, 2, 2))
  boxplot(data)
  if(file != ""){
    dev.off()
  }
}

########main procedure########

setwd("H:/shengquanhu/projects/chenxi/20130724-chenxi_colon_cancer");
genemap<-read.table("H:/shengquanhu/projects/database/affy/HG-U133_Plus_2.probe_gene_map.txt", header=T, as.is=T, sep="\t", row.names=1)

training<-"coloncancer_training.csv"
stage23<-"coloncancer_training_stage23.csv"
if(!file.exists(stage23)){
  cat(paste0("reading file ", training, " ...\n"))
  
  #read/transfer data
  trainingdata<-read.csv(training, header=T, row.names=1, check.names=F)
  
  sample<-read.table("data/data_SampleInformation.txt", sep="\t", header=T, row.names=2)
  
  isstage2<-str_detect(sample$Stage, "^[2]") | sample$Stage == "AJCC stage II CRC"
  isstage3<- str_detect(sample$Stage, "^[3]")
  
  isstage2<-isstage2 | ((!isstage3) & sample$DukeStage == 'B')
  isstage3<-isstage3 | ((!isstage2) & sample$DukeStage == 'C')
  
  isstage23<-isstage2 | isstage3
  
  table(isstage23)
  
  sample$FinalStage<-""
  sample$FinalStage[isstage2]<-"2"
  sample$FinalStage[isstage3]<-"3"
  
  stage23sample<-sample[isstage23,]
  
  instage23<-rownames(trainingdata) %in% rownames(stage23sample)
  stage23training<-trainingdata[instage23, ]
  
  cat(paste0("writing file ", stage23, " ...\n"))
  write.csv(stage23training, file=stage23)
  
  batchdefinitionfile<-paste0(training,".batchdefinition");
  batchdata<-read.csv(batchdefinitionfile, header=T)
  
  stage23batchdata<-batchdata[batchdata$Sample %in% rownames(stage23training),]
  stage23batchdata$Stage<-stage23sample[as.character(stage23batchdata$Sample),]$FinalStage
  
  stage23batchdefinitionfile<-paste0(stage23,".batchdefinition");
  write.csv(stage23batchdata, file=stage23batchdefinitionfile,row.names=F)
}

datafiles<-c(
  #"coloncancer_cellline.tsv",
  #"coloncancer_testing.tsv"
  stage23,
  training
)

datafile<-datafiles[1]

for(datafile in datafiles){
  #quantile normalization
  qnormfile<-paste0(datafile,".qnorm")
  qnormfileR<-paste0(qnormfile,".RData")
  if(!file.exists(qnormfileR)){
    cat(paste0("reading file ", datafile, " ...\n"))
    
    #read/transfer data
    data<-read.csv(datafile, header=T, row.names=1, check.names=F)
    tdata<-log2(t(data))
    rm(data)
    
    drawBoxplot(tdata, file=paste0(datafile, ".boxplot.png"), width=20000, height=2000)
    
    cat("quantile normalization ...\n")
    
    normalize.quantiles.robust(x=tdata,copy=F)
    save(tdata,file=qnormfileR)
    
    drawBoxplot(tdata, file=paste0(qnormfile, ".boxplot.png"), width=20000, height=2000)
  }
  
  parseGenes(qnormfile)
  
  bnormfile<-paste0(qnormfile,".bnorm")
  bnormfileR<-paste0(bnormfile,".RData")
  if(!file.exists(bnormfileR)){
    load(qnormfileR)
    
    cat("batch normalization by ComBat function ...\n")
    
    batchdefinitionfile<-paste0(datafile,".batchdefinition");
    batchdata<-read.csv(batchdefinitionfile, header=T)
    batch<-batchdata$Batch
    
    if(length(unique(batch)) > 1){
      mod<-model.matrix(~1, data=batchdata)
      bdata<-ComBat(dat=tdata, batch=batch, mod=mod)
    }else{
      bdata<-tdata
    }      
    rm(tdata)
    save(bdata,file=bnormfileR);
    
    drawBoxplot(bdata, file=paste0(bnormfile, ".boxplot.png"), width=20000, height=2000)
  }else{
    load(bnormfileR)
  }
  
  ret<-parseGenes(bnormfile)
  
  require("WGCNA")
  options(stringsAsFactors = FALSE);
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  
  load(ret[4])
  datExpr<-t(pregdata)
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  png(paste0(ret[3], ".findpower.png"), width=4000, height=3000,res=300)
  
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  dev.off()
  
  #based on the image, 6 is the power
  powertolerance=6
  
  tpregdata<-t(pregdata)
  net = blockwiseModules(tpregdata, power = powertolerance,
                         TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "networkPreg",
                         verbose = 3)
  
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  png(paste0(ret[3], ".cluster.png"), width=4000, height=3000,res=300)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  
  save(MEs, moduleLabels, moduleColors, geneTree, file = paste0(ret[3], ".netauto.RData"))
  
  # Calculate topological overlap anew: this could be done more efficiently by saving the TOM
  # calculated during module detection, but let us do it again here.
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = powertolerance);
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM = dissTOM^7;
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA;
  # Call the plot function
  png(paste0(ret[3], ".heatmap.png"), width=4000, height=3000,res=300)
  TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  dev.off()
}

cat("done\n")
