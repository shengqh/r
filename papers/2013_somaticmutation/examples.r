library("brglm") 
library("elrm")

file<-"H:/shengquanhu/projects/somaticmutation/rsmc/result/TCGA-A7-A0D9/candidates/StrandFailed_6_31323974_C.wsm"
file<-"H:/shengquanhu/projects/vangard/VANGARD00095_jennifer_exome/rsmc/StrandFailed_1_21809750_T.wsm" #score sample,
file<-"H:/shengquanhu/projects/vangard/VANGARD00095_jennifer_exome/rsmc/PositionFailed_5_180041208_T.wsm" #position sample,


#get chr and position from filename
parts<-unlist(strsplit(file,'_'))
chr<-parts[2]
position<-parts[3]
parts<-unlist(strsplit(parts[4],'[.]'))
nucleotide<-parts[1]

raw.data<-read.table(file, header=T)
tb<-table(raw.data$Base)
if(tb[1] < tb[2]){
  levels=names(tb)
}else{
  levels=c(names(tb)[2], names(tb)[1])
}

#  fit<-brglm(Base ~ factor(SAMPLE) + Score + Strand, family=binomial, data=raw.data) 
#fit<-brglm(Base ~ factor(SAMPLE) + Strand, family=binomial, data=raw.data) 
fit<-brglm(Base ~ factor(SAMPLE) + Score + Position, family=binomial, data=raw.data) 
fit.coef=coef(summary(fit))
fit.coef
