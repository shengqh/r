library("VennDiagram")

setwd("H:/shengquanhu/projects/somaticmutation/")

getintersect<-function(data, names){
  result<-data
  lapply(names, function(name){
    result<<-result[result[,name] != '',]
  })
  nrow(result)
}

drawvenn<-function(data, cnames){
  plot.new()
  draw.quad.venn(getintersect(data, cnames[1]),
                 getintersect(data, cnames[2]),
                 getintersect(data, cnames[3]),
                 getintersect(data, cnames[4]),
                 
                 getintersect(data, cnames[c(1:2)]),
                 getintersect(data, cnames[c(1,3)]),
                 getintersect(data, cnames[c(1,4)]),
                 getintersect(data, cnames[c(2,3)]),
                 getintersect(data, cnames[c(2,4)]),
                 getintersect(data, cnames[c(3,4)]),

                 getintersect(data, cnames[c(1,2,3)]),
                 getintersect(data, cnames[c(1,2,4)]),
                 getintersect(data, cnames[c(1,3,4)]),
                 getintersect(data, cnames[c(2,3,4)]),
                 
                 getintersect(data, cnames[c(1,2,3,4)]),
                 
                 category=cnames
  )
}

calculateTable<-function(data, cnames){
  refcount<-getintersect(data, cnames[1])
  res<-matrix(unlist(lapply(c(2:length(cnames)), function(i){
    methodcount<-getintersect(data, cnames[i])
    overlapcount<-getintersect(data, cnames[c(1,i)])
    c(refcount, methodcount, overlapcount, overlapcount / refcount, overlapcount / methodcount)
  })), ncol=5, byrow=TRUE)
  colnames(res)<-c(cnames[1], "Count", "OverlapCount", "Sensitivity", "Specificity")
  rownames(res)<-c(cnames[2:length(cnames)])
  res
}

all<-read.table("TCGA_muTect_varscan2_rsmc.site.tsv",sep='\t',header=T)

pdf("TCGA_muTect_varscan2_rsmc.site.DNA.pdf")
drawvenn(all, colnames(all)[c(4:7)])
dev.off()

pdf("TCGA_muTect_varscan2_rsmc.site.RNA.pdf")
drawvenn(all, colnames(all)[c(4,8:10)])
dev.off()

read10<-read.table("TCGA_muTect_varscan2_rsmc.site.depth10",sep='\t',header=T)
matched<-all[read10$IsDepth10 == "True",]

matched$DNA_M_V<-(matched$DNA_MUTECT != '' & matched$DNA_VARSCAN2 != '')
matched$DNA_M_V[!matched$DNA_M_V]<-''

pdf("TCGA_muTect_varscan2_rsmc.site.depth10.DNA.pdf")
drawvenn(matched, colnames(matched)[c(4:7)])
dev.off()

pdf("TCGA_muTect_varscan2_rsmc.site.depth10.RNA.pdf")
drawvenn(matched, colnames(matched)[c(4,8:10)])
dev.off()

write.csv(calculateTable(matched, colnames(matched)[c(5:8)]), file="TCGA_muTect_varscan2_rsmc.site.DNA.depth10.csv")

write.csv(calculateTable(matched, colnames(matched)[c(5,9:11)]), file="TCGA_muTect_varscan2_rsmc.site.RNA.depth10.csv")

pdf("H:/shengquanhu/projects/somaticmutation/TCGA_muTect_varscan2_rsmc.site.RNA_MV.depth10.pdf")
drawvenn(matched, colnames(matched)[c(19,9:11)])
dev.off()

write.csv(calculateTable(matched, colnames(matched)[c(19,9:11)]), file="TCGA_muTect_varscan2_rsmc.site.RNA_MV.depth10.csv")

