library("affy");
library("sva")

setwd("H:/shengquanhu/projects/chenxi/20130724-chenxi_colon_cancer");

sampleInfo<-read.csv("sampleinfo.csv", header=T, check.names=F, row.names=1)
stage23Info<-sampleInfo[!is.na(sampleInfo$FinalStage),]
rm(sampleInfo)

stage23Names<-unique(as.character(stage23Info$Dataset))

sampleWithSurvival<-stage23Info[!is.na(stage23Info$DFSurvivalTime),]
survivalNames<-unique(as.character(sampleWithSurvival$Dataset))
survivalNames<-survivalNames[! (survivalNames %in% "GSE39582")]

datasets<-list()
datasets[["training"]]<-c(as.character(stage23Names[!(stage23Names %in% survivalNames)]))
for(geo in survivalNames){
  datasets[[geo]]<-c(geo)
}

dataTable<-data.frame(GeoName = c(datasets[["training"]], survivalNames[survivalNames != "GSE39582"]),
                      GeoType = c("survival/training", rep("training", length(datasets[["training"]]) - 1), rep("testing", length(datasets)-2)))

dataTable$GsmCount<-unlist(lapply(dataTable$GeoName, function(x){
  nrow(stage23Info[stage23Info$Dataset %in% x,])
}))

write.csv(dataTable, "stage23table.csv", row.names=F)

#dsName<-"training"
for(dsName in names(datasets)){
  print(dsName)
  ds<-datasets[[dsName]]
  dataInfo<-stage23Info[stage23Info$Dataset %in% ds,]
  
  rmaFilename<-paste0(dsName, ".rma.RData")

  if(exists("rmaData")){
    rm(rmaData)
  }
  
  if(!file.exists(rmaFilename)){
    sampleFiles<-as.character(dataInfo$SampleFile)
    sampleFiles<-unlist(lapply(sampleFiles,function(x){
      substring(x, nchar(getwd())+1)
    }))
    rmaData<-justRMA(filenames=sampleFiles, sampleNames=rownames(dataInfo))
    rmaExpr<-exprs(rmaData)
    save(rmaExpr, dataInfo, file=rmaFilename)
    rm(rmaData)
  }
  
  batch<-dataInfo$BatchEffect
  if(length(unique(batch)) > 1){
    rmaCombatFilename<-paste0(dsName, ".rma.combat.RData")
    if(!file.exists(rmaCombatFilename)){
      if(!exists("rmaExpr")){
        load(rmaFilename)
      }
      
      mod<-model.matrix(~1, data=dataInfo)
      
      batchExpr<-ComBat(dat=rmaExpr, batch=batch, mod=mod)
      save(batchExpr, dataInfo, file=rmaCombatFilename)
      rm(rmaExpr)
      rm(batchExpr)
    }
  }
}
