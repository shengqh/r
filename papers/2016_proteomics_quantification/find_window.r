setwd("H:/shengquanhu/projects/20151103_shilin_labelfree_elution_order_algorithm/test")

loadOrInstallPackage <- function(x)
{
  if (!require(x,character.only = TRUE, quietly = TRUE))
  {
    options("repos"="http://cran.us.r-project.org")
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE, quietly = TRUE)) {
      stop(paste0("Package ", x, " not found"))
    }
  }
}

loadOrInstallPackage("ggplot2")
loadOrInstallPackage("changepoint")
loadOrInstallPackage("cowplot")

MAX_ISOTOPIC<-3
MINIMUM_MEAN_PPM_DIFF<-2
MINIMUM_VARIANCE_PPM_DIFF<-2
MINIMUM_SCAN_EACH_SIDE<-5

findBoundary<-function(isodata, maxIdentifiedIndex, cpt, minDiff, minScan){
  leftBoundary<-1
  if(maxIdentifiedIndex > minScan){
    res<-cpt(isodata$PPMTolerance[1:maxIdentifiedIndex])
    if(ncpts(res) > 0){
      diffmean<-abs(param.est(res)$mean)
      resmean<-param.est(res)$mean
      if(abs(resmean[1]) > abs(resmean[2])){
        diffmean<-abs(resmean[1] -resmean[2])
        if(diffmean > minDiff){
          leftBoundary<-max(1, min(maxIdentifiedIndex-minScan, cpts(res)))
        }
      }
    }
  }
  
  rightBoundary<-nrow(isodata)
  if(nrow(isodata) - maxIdentifiedIndex > minScan){
    res<-cpt(isodata$PPMTolerance[maxIdentifiedIndex:nrow(isodata)])
    if(ncpts(res) > 0){
      diffmean<-abs(param.est(res)$mean)
      resmean<-param.est(res)$mean
      if(abs(resmean[1]) < abs(resmean[2])){
        diffmean<-abs(resmean[1] -resmean[2])
        if(diffmean > minDiff){
          rightBoundary<-min(nrow(isodata), max(maxIdentifiedIndex + cpts(res), maxIdentifiedIndex + minScan))
        }
      }
    }
  }
  
  boundaries<-isodata$RetentionTime[c(leftBoundary, rightBoundary)]
  return (boundaries)
}

df <- data.frame(File=character(), 
                 LeftBoundary=numeric(), 
                 RightBoundary=numeric(),
                 stringsAsFactors = FALSE)
#inputFile="CT102_0_Hour_DLGEENFK_476_10160.chro.tsv"
#inputFile="CT102_0_Hour_LGQYASPTAK_518_6738.chro.tsv"
#inputFile="CT102_0_Hour_LSQRFPK_438_5534.chro.tsv"
#inputFile="CT102_0_Hour_LSPLGEEMRDR_435_10897.chro.tsv"
inputFile="CT102_0_Hour_QTALVELVK_501_16906.chro.tsv"

cat(inputFile, "/n")
out<-read.table(inputFile, header=T)
out$Isotopic = factor(out$Isotopic)

iso0data<-out[out$Isotopic == "1",]
identified<-iso0data[iso0data$Identified == "True",]
if(nrow(identified) == 0){
  maxIdentifiedIndex=round(nrow(iso0data) / 2)
}else{
  index <-which.max(identified$Intensity)
  maxIdentifiedIndex<-which(iso0data$Scan==identified$Scan[index])
}
maxIdentifiedRetentionTime<-iso0data$RetentionTime[maxIdentifiedIndex]

isos<-unique(out$Isotopic)
isoBoundary<-out$RetentionTime[c(1,nrow(out))]
for(iso in isos[1:MAX_ISOTOPIC]){
  isodata<-out[out$Isotopic == iso,]
  curBoundary<-findBoundary(isodata, maxIdentifiedIndex, cpt.mean, MINIMUM_MEAN_PPM_DIFF, MINIMUM_SCAN_EACH_SIDE)
  isoBoundary[1] = max(isoBoundary[1], curBoundary[1])
  isoBoundary[2] = min(isoBoundary[2], curBoundary[2])
  curBoundary<-findBoundary(isodata, maxIdentifiedIndex, cpt.meanvar, MINIMUM_VARIANCE_PPM_DIFF, MINIMUM_SCAN_EACH_SIDE)
  isoBoundary[1] = max(isoBoundary[1], curBoundary[1])
  isoBoundary[2] = min(isoBoundary[2], curBoundary[2])
}

newdata<-data.frame(File=inputFile, 
                    LeftBoundary=isoBoundary[1], 
                    RightBoundary=isoBoundary[2])
df<-rbind(df, newdata)

iso0data<-out[out$Isotopic == 1,c("RetentionTime", "Scan", "Intensity", "Identified")]
identified<-iso0data[iso0data$Identified == "True",]
if(nrow(identified) == 0){
  maxIdentifiedIndex=round(nrow(iso0data) / 2)
}else{
  index <-which.max(identified$Intensity)
  maxIdentifiedIndex<-which(iso0data$Scan==identified$Scan[index])
}
maxIdentifiedRetentionTime<-iso0data$RetentionTime[maxIdentifiedIndex]

identifiedScan = unique(out$RetentionTime[out$Identified == "True"])

iso0win<-iso0data[iso0data$RetentionTime >= isoBoundary[1] & iso0data$RetentionTime <= isoBoundary[2],]
maxIndex<-which.max(iso0win$Intensity)
maxRetention=iso0win$RetentionTime[maxIndex]
minIntensity<-iso0win$Intensity[maxIndex] * 0.1

removeLeft<-iso0win[iso0win$Intensity < minIntensity & iso0win$RetentionTime < maxRetention,]
removeRight<-iso0win[iso0win$Intensity < minIntensity & iso0win$RetentionTime > maxRetention,]

newwinLeft<-ifelse(nrow(removeLeft) == 0, iso0win$RetentionTime[1], max(removeLeft$RetentionTime))
newwinRight<-ifelse(nrow(removeRight) == 0, iso0win$RetentionTime[nrow(iso0win)], min(removeRight$RetentionTime))

#out<-out[out$RetentionTime > 54,]

p2<-ggplot(out, aes(x=RetentionTime, y=Intensity, group=Isotopic)) + geom_line(aes(color=Isotopic))
if(nrow(identified) != 0){
  p2<-p2+
    geom_vline(xintercept = identifiedScan, colour="black") +
    geom_vline(xintercept = maxIdentifiedRetentionTime, colour="red")
}
p2<-p2+
  geom_vline(xintercept = isoBoundary, colour="blue") +
  #geom_vline(xintercept = c(newwinLeft, newwinRight), colour="green") +
  theme_cowplot()
print(p2)

