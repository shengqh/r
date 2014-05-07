library("stringr")

setwd("H:/shengquanhu/projects/chenxi/20130724-chenxi_colon_cancer");

tissueinfo<-read.table("data/data_SampleInformation.tsv", header=T, check.names=F, row.names=2, sep="\t")

isstage2<-str_detect(tissueinfo$Stage, "^[2]") | tissueinfo$Stage == "AJCC stage II CRC"
isstage3<- str_detect(tissueinfo$Stage, "^[3]")

isstage2<-isstage2 | ((!isstage3) & tissueinfo$DukeStage == 'B')
isstage3<-isstage3 | ((!isstage2) & tissueinfo$DukeStage == 'C')

isstage23<-isstage2 | isstage3

tissueinfo$FinalStage<-""
tissueinfo$FinalStage[isstage2]<-"2"
tissueinfo$FinalStage[isstage3]<-"3"

tissueinfo$BatchEffect<-unlist(apply(tissueinfo, 1, function(x){
  paste0(x["Dataset"], "_", x["Batch"])
}))

write.csv(tissueinfo, "sampleinfo.csv")
