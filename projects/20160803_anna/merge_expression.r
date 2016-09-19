setwd("H:/shengquanhu/projects/20160803_anna_proteomics_metabolomics")

d0<-read.delim("proteomics_expression_data.txt", header=T, stringsAsFactors = F)
d0<-d0[,c(3,12:ncol(d0))]
d0$ID_<-gsub("\\s+","", d0$ID_)
d0$ID_<-gsub("C15T0061212","C15-T0061212", d0$ID_)

info<-read.csv("Brittain Metadata RTI.csv")
rownames(info)=info$ID
sampleinfo=info[d0$ID_,]
sampleinfo["PPH189CM5264", "PAH"] = 1
sampleinfo$Group="Unknown"
sampleinfo$Group[sampleinfo$NoPH==1] = "NoPH"
sampleinfo$Group[sampleinfo$PVH==1] = "PVH"
sampleinfo$Group[sampleinfo$Combined==1] = "Combined"
sampleinfo$Group[sampleinfo$PAH==1] = "PAH"
sampleinfo$Group[sampleinfo$CTD.PAH==1] = "CTD-PAH"
sampleinfo$Group[sampleinfo$UMC==1] = "UMC"
sampleinfo$Group=as.factor(sampleinfo$Group)
sampleinfo=sampleinfo[,c("ID", "Group")]

table(sampleinfo$Group)

write.table(sampleinfo, file="sample.group", row.names=F, sep="\t", quote=F)

library(heatmap3)

hmcols <- colorRampPalette(c("green", "black", "red"))(256)
drawHCA<-function(prefix, rldselect, conditionColors, groupcols, gnames){
  htfile<-paste0(prefix, "-heatmap.png")
  cat("saving HCA to ", htfile, "\n")
  genecount<-nrow(rldselect)
  if(genecount > 2){
    png(filename=htfile, width=3000, height=3000, res=300)
    heatmap3(rldselect,
             col = hmcols,
             ColSideColors = conditionColors,
             margins=c(12,5),
             scale="r",
             dist=dist,
             labRow=NA,
             main=paste0("Hierarchical Cluster Using ", genecount, " Genes"),
             useRaster=FALSE,
             legendfun=function() showLegend(legend=gnames, col=groupcols,cex=1.0,x="center"))
    dev.off()
  }
}

files=c("proteomics_multivariable_PAH_vs_Combined_winners",
        "proteomics_multivariable_PAH_vs_NoPH_winners",
        "proteomics_multivariable_PVH_vs_Combined_winners",
        "proteomics_multivariable_PVH_vs_PAH_winners")

file=files[1]
allwinner=c()
for(file in files){
  winner=read.csv(paste0(file, ".csv"), header=T, stringsAsFactors = F)
  winnerdata=d0[, c("ID_", winner$X)]
  colnames(winnerdata)[1] = "Sample"
  write.csv(file=paste0(file,".expression.csv"), winnerdata, row.names=F)
  allwinner=c(allwinner, winner$X)
}

unique_allwinner=unique(allwinner)
allwinnerdata=d0[, unique_allwinner]
rownames(allwinnerdata)=d0$ID_
allwinnerdata=t(allwinnerdata)

countNum=log(allwinnerdata)

groupcols=c("red", "blue", "green", "yellow" )
conditionColors<-as.matrix(data.frame(Group=groupcols[sampleinfo$Group]))
drawHCA("winner", countNum, conditionColors = conditionColors, groupcols, levels(sampleinfo$Group) )

PAH_NoPH_winner=read.csv(paste0(files[2], ".csv"), header=T, stringsAsFactors = F)
subsample=sampleinfo[sampleinfo$Group=="PAH" | sampleinfo$Group=="NoPH",]
subdata=countNum[PAH_NoPH_winner$X, ]

subdata=countNum[PAH_NoPH_winner$X, as.character(subsample$ID)]
subcolors<-as.matrix(data.frame(Group=groupcols[subsample$Group]))
subgroups=c("NoPH", "PAH")
subgroupcolors=c("blue", "green")
conditionColors<-as.matrix(data.frame(Group=groupcols[subsample$Group]))
drawHCA("PAH_NoPH_winner", subdata, conditionColors = conditionColors, groupcols=subgroupcolors, gnames=subgroups  )

