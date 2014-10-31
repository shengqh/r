ppm2mz<-function(precursorMz, ppmTolerance)
{
  precursorMz * ppmTolerance / 1000000;
}

mz2ppm<-function(precursorMz, mzTolerance)
{
  mzTolerance * 1000000 / precursorMz;
}

setwd("E:/shengquanhu/projects/suzhiduan")
data<-read.table("20140829_Q1_Tim_ROS_iodoTMT6_F1.raw.oitraq.xml.param", sep="\t", header=T)

probabilities<-apply(data, 1, function(x){
  mzt<-ppm2mz(x["Ion"], 10)
  pnorm(x["Ion"]+mzt, x["DetectedIon"], x["DetectedSigma"]) - pnorm(x["Ion"]-mzt, x["DetectedIon"], x["DetectedSigma"])
})

data$IonWithin10ppm<-probabilities
