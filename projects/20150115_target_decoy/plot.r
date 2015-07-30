library(ggplot2)

setwd("E:/shengquanhu/projects/20150115_percolator/target_decoy_spectra/peptides")

datasets<-c("Fusion_HCDIT_Yeast","Fusion_HCDOT_Human","QExactive_HCDOT_Human","QTOF_Ecoli", "Elite_CIDIT_Human", "Fusion_CIDIT_Human")
dataset<-datasets[1]

alldata<-NULL
for(dataset in datasets){
  cat(dataset, "\n")
  
  pngfile<-paste0(dataset,"_MSGF.expectvalue.png")
  
  #if(!file.exists(pngfile)){
    data1<-read.table(paste0(dataset, ".plus0.1dalton.msgf.mzid.peptides"), sep='\t', header=T)
    data1$Scan<-as.numeric(sub("^[^,]*,", "", data1$FileScan))
    data1$Type<- ifelse(data1$Scan > 3000000, "Shifted 0.1 Dalton", "Unshifted")
    
    data2<-read.table(paste0(dataset,".plus10dalton.msgf.mzid.peptides"), sep='\t', header=T)
    data2$Scan<-as.numeric(sub("^[^,]*,", "", data2$FileScan))
    data2<-data2[data2$Scan > 3000000,]
    data2$Type<-"Shifted 10 Daltons"
    
    data<-rbind(data1, data2)
    
    if(is.null(alldata)){
      alldata<-data
      alldata$Dataset<-dataset
    }else{
      curdata<-data
      curdata$Dataset<-dataset
      alldata<-rbind(alldata, curdata)
    }
    
    data$logEvalue<--log(data$ExpectValue)
    data$Database<-ifelse(grepl("REVERSE_",data$Reference), "Decoy", "Target")
    
    data<-data[data$Charge < 5,]
    data$Charge<-paste0("Charge ", data$Charge, "+")
    
    png(pngfile, width=3000, height=3000, res=300)
    g<-ggplot(data, aes(Database, logEvalue)) + geom_violin(aes(fill=Database)) + facet_grid(Type ~ Charge) + ylab("-log(ExpectValue)") + ggtitle(dataset)
    print(g)
    dev.off()
  #}
}

pie<-ggplot(alldata, aes(x=factor(1), fill=Charge)) + 
  geom_bar(width=1) + 
  coord_polar(theta="y") + 
  facet_grid(. ~ Dataset)
print(pie)