setwd("H:/shengquanhu/projects/20160803_anna_proteomics_metabolomics")

files=c("proteomics_multivariable_PAH_vs_Combined_winners",
        "proteomics_multivariable_PAH_vs_NoPH_winners",
        "proteomics_multivariable_PVH_vs_Combined_winners",
        "proteomics_multivariable_PVH_vs_PAH_winners")

d0<-read.delim("proteomics_expression_data.txt", header=T, stringsAsFactors = F)
d0<-d0[,c(3,12:ncol(d0))]
d0$ID_<-gsub("\\s+","", d0$ID_)
d0$ID_<-gsub("C15T0061212","C15-T0061212", d0$ID_)
d1<-t(d0)

file=files[1]
for(file in files){
  winner=read.csv(paste0(file, ".csv"))
  winnerdata=d1[winner$X,]
  merged=cbind(winner, winnerdata)
  colnames(merged)=c("protein", "multivariable_value", d1["ID_",])
  write.csv(file=paste0(file,".expression.csv"), merged, row.names=F)
}