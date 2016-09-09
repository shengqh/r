setwd("H:/shengquanhu/projects/20160830_liuqi_genepanel/documents")
data<-read.csv("TCPS_Labcorp_Link_20160829.csv", stringsAsFactors = F)

data$tmp<-data$recur_status
data$tmp[data$tmp == "No recur"] = "Norecur"
data$ID = paste(data$StudyID, data$tmp, data$BaseRecur, sep="_")

data<-data[,c("barcode", "ID", "StudyID", "BaseRecur", "status","recur_status")]
write.table(data, file="TCPS_Labcorp_Link_20160829.txt", row.names = F,col.names =T,  quote = F, sep="\t")
