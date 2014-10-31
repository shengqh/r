
library("outliers")

findOutlier<-function(inputFile, outputFile, recursive=TRUE, confidence=0.99, minFinalCount=3, subjectIndex=1, ratioIndex=2){
  pvalue<-1-confidence
  
  data<-read.table(inputFile, header=T, sep="\t")
  
  data$IsOutlier<-rep(FALSE, nrow(data))
  
  subjects<-unique(data[,subjectIndex])
  
  result<-data[FALSE,]
  for(sub in subjects){
    subdata<-data[data[,subjectIndex] == sub,]  
    while(nrow(subdata) > minFinalCount){
      logratio<-log(subdata[,ratioIndex])
      gt<-grubbs.test(logratio)
      if(gt$p.value < pvalue ){
        if(grepl("lowest",gt$alternative )){
          outlierindex <- which.min(logratio)
        }else{
          outlierindex<-which.max(logratio)
        }
        subdata[outlierindex,"IsOutlier"]<-TRUE
        result<-rbind(result,subdata[outlierindex,])
        subdata<-subdata[!subdata$IsOutlier,]
        
        if(!recursive){
          break
        }
      }else{
        break
      }
    }
    result<-rbind(result, subdata)
  }
  
  result<-result[rownames(data),]
  
  write.table(result, outputFile, row.names=F, sep="\t")
}

setwd("E:/sqh/Dropbox/career/papers/draft/2014-iCan")
findOutlier("MAP_frames_M3P001.txt", "outlier.tsv")


