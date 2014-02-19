library("cn.mops")

cnmops<-load("H:/shengquanhu/projects/CNV/cnmops/result/2110_resCNMOPS_exomecn.mops.Rdata")
cnvs<-resCNMOPS@cnvs

d<-cbind(substring(as.character(cnvs@elementMetadata@listData$sampleName),2),
         as.character(cnvs@seqnames), 
         as.character(cnvs@ranges@start),
         as.character(as.numeric(cnvs@ranges@start) + as.numeric(cnvs@ranges@width) - 1),
         as.character(cnvs@ranges@width),
         as.character(cnvs@elementMetadata@listData$CN),
         as.character(cnvs@elementMetadata@listData$median),
         as.character(cnvs@elementMetadata@listData$mean))
colnames(d)<-c("sample","chr","start","end", "length","type","median","mean")
d<-d[order(d[,"sample"], as.numeric(d[,"chr"]), as.numeric(d[,"start"])),]
d[,"chr"]<-paste0("chr",d[,"chr"])
d[,"type"]<-apply(d,1,function(x){
  if(as.numeric(x["median"]) < 0){
    return ("DELETION")
  }else{
    return ("DUPLICATION")
  }
})

write.table(d, file="H:/shengquanhu/projects/CNV/cnmops.txt",sep="\t",col.names=T,row.names=F,quote=F)
