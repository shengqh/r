setwd("H:/shengquanhu/projects/Jennifer/20150805_johanna_neuroblastoma")

library(limma)  # two-color pre-processing; differential
library(hgu133plus2.db)
library(heatmap3)

evaluesfile<-"GSE49710.rdata"
if(file.exists(evaluesfile)){
  load(evaluesfile)
}else{
  targets <- readTargets("GSE49710/targets.txt")
  
  x <- read.maimages(targets, path="GSE49710", source="agilent",green.only=TRUE)
  y <- backgroundCorrect(x, method="normexp", offset=16)
  y <- normalizeBetweenArrays(y, method="quantile")
  
  save(y, file=evaluesfile)
}
