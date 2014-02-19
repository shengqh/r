setwd("E:/sqh/Dropbox/papers/published/2013-CNV/paper")

library(ggplot2)
library(reshape2)
library(flexmix)

readdata<-function(filename, method){
  data<-read.table(filename,sep="\t", header=T)
  tt<-table(data$sample, data$type)
  mt<-melt(tt,id.vars=1:2)
  mtt<-cbind(c(rep(method,nrow(mt))), mt)
  return(mtt)
}

names<-c("aCGH","cn.MOPS","CoNIFER","exomeCopy", "ExomeDepth")
filenames<-c("cgh.txt","cnmops.txt","conifer.txt","exomeCopy.txt","ExomeDepth.txt")
ncount<-length(names)

deletion<-c(1:16)
duplication<-c(17:32)

data<-lapply(c(1:ncount), FUN=function(x){
  return (readdata(filenames[x],names[x]))
})

mt2<-cbind(data[[1]]$value[deletion], 
           data[[2]]$value[deletion], 
           data[[3]]$value[deletion], 
           data[[4]]$value[deletion], 
           data[[5]]$value[deletion]
)
colnames(mt2)<-names
kd<-KLdiv(mt2, method="discrete")
rk<-c(1,2,4,5,3)
kd<-kd[rk,rk]
write.csv(format(kd, digits=2), file="deletion.similarity.csv")

mt2<-cbind(data[[1]]$value[duplication], 
           data[[2]]$value[duplication], 
           data[[3]]$value[duplication], 
           data[[4]]$value[duplication], 
           data[[5]]$value[duplication])
colnames(mt2)<-names
kd<-KLdiv(mt2, method="discrete")
kd<-kd[rk,rk]
write.csv(format(kd, digits=2), file="duplication.similarity.csv")

pvalues<-lapply(c(1:ncount), FUN=function(x){
  return (wilcox.test(x=data[[x]]$value[deletion],y=data[[x]]$value[duplication],paired=TRUE)$p.value)
})

pvalues.adjust<-p.adjust(pvalues, method="BH")
pvaluestrs<-paste0("p=",formatC(pvalues.adjust, digits=2))

data<-lapply(c(1:ncount), FUN=function(x){
  return (readdata(filenames[x],paste0(names[x], " [", pvaluestrs[x], "]")))
})

mtt<-rbind(data[[1]], data[[4]],data[[2]], data[[5]], data[[3]])
colnames(mtt)<-c("Method","Sample","CNV","Count")

pdf(file="figure2.pdf", width=6, height=9)

ggplot(mtt, aes(x=Sample, y=Count, fill=CNV)) + 
  geom_bar(stat="identity", position="stack") + 
  facet_wrap(~Method, ncol=1, scales="free_y") + theme_bw() + 
  theme(axis.text.x=element_text(angle=90, vjust = 0.5))

gap<-0.174
y<-0.963-gap*c(0:4)
grid.text(label = c("(A)","(B)","(C)","(D)","(E)"), x = 0.005, y = y, just = c(0,0), gp = gpar(col = "black", cex = 1.5, fontface = "plain", fontfamily = "serif"))

dev.off()

# data2<-lapply(c(1:ncount), FUN=function(x){
#   tf<-lapply(c(1:16), FUN=function(y){
#     if(data[[x]][y+16,4] > data[[x]][y,4]){
#       return (data[[x]][y+16,3])
#     }else{
#       return (data[[x]][y,3])
#     }
#   })
#   return (cbind(names[x], t(table(unlist(tf)))))
# })
# dup_del<-rbind(data2[[1]][,c(1,3,2)], 
#                data2[[2]][,c(1,3,2)], 
#                data2[[3]][,c(1,3,2)], 
#                data2[[4]][,c(1,3,2)], 
#                data2[[5]][,c(1,3,2)],
# )
# colnames(dup_del)[1]<-"Method"
# write.table(dup_del, file="Duplication_Deletion_Sample.tsv",sep="\t",quote=F,row.names=F)
