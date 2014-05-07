


library(hgu133plus2.db)
library(affy)
library(annotate)
library(survival)
library(multtest)
library(GEOquery)

## aa<-getGEO("GSE39582")
## GSE39582<-aa
## save(GSE39582, file="GSE39582.RData")
load("GSE39582.RData")
pinfo<-pData(phenoData(GSE39582[[1]]))

pinfo<-pinfo[,c(11, 12, 13,16,17, 18, 24,28,32)]
colnames(pinfo)<-c("sex", "age", "stage", "rfs_event", "rfs_surtime", "MMR", "BRAF", "KRAS", "TP53")
pinfo$sex<-gsub("Sex: ", "", pinfo$sex)
pinfo$age<-gsub("age.at.diagnosis: ", "", pinfo$age)

pinfo$stage<-gsub("tnm.stage: ", "", pinfo$stage)
pinfo$rfs_event<-gsub("rfs.event: ", "", pinfo$rfs_event)
pinfo$rfs_surtime<-gsub("rfs.delay: ", "", pinfo$rfs_surtime)
pinfo$MMR<-gsub("mmr.status: ", "", pinfo$MMR)
pinfo$BRAF<-gsub("braf.mutation.status: ", "", pinfo$BRAF)
pinfo$KRAS<-gsub("kras.mutation.status: ", "", pinfo$KRAS)
pinfo$TP53<-gsub("tp53.mutation.status: ", "", pinfo$TP53)

array<-exprs(GSE39582[[1]])


probelist<-rownames(array)

geneSymbol<-getSYMBOL(probelist,"hgu133plus2.db")

anno<-data.frame(probelist, geneSymbol)

pinfo1<-pinfo[(pinfo$stage==2 | pinfo$stage==3),]
pinfo1<-pinfo1[(pinfo1$rfs_event!="NA" & pinfo1$rfs_surtime!="NA"),]
pinfo1$age<-as.numeric(pinfo1$age)
array1<-array[, rownames(pinfo1)]

### survival analysis


survF=function(gene.expression,respSurv){ 
  tmp=summary(coxph(respSurv~gene.expression))
  return(c(tmp$coef[1,])
           )

}

## OS
objSurv=Surv(as.numeric(pinfo1$rfs_surtime), as.numeric(pinfo1$rfs_event))
res=t(apply(array1,1,survF,respSurv=objSurv)) 
dim(res)



colnames(res)=c("coef", "exp_coef", "se", "z", "pvalue")
head(res)

bh<-mt.rawp2adjp(res[,5], "BH")
adjustedP<-bh$adjp[order(bh$index),]

res<-data.frame(res, adjustedP)
res<-res[,-6]
res<-data.frame(res, anno$geneSymbol)

write.csv(res, "stage23.csv")

