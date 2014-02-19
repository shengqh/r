source("e:/sqh/programs/r/mirna/chris_v2/common.r")

brca<-loadData("TCGA_mirnaseq_rnaseqv2_20130521_brca_rank_sample.csv")

colnames(brca)[1:5]<-c("mir141", "mir200a", "mir200b", "mir200c", "mir429")

mirna<-brca[,c(1:5)]
rnaseqv2<-brca[,c(6:ncol(brca))]

glmres<-apply(rnaseqv2, 2, function(x){
  summ<-summary(glm(x ~ mirna[,1] + mirna[,2] + mirna[,3] + mirna[,4] + mirna[,5]))
  coeff<-coefficients(summ)
  return (coeff[,4])
})

glmres<-t(glmres)
saveData(glmres, "TCGA_mirnaseq_rnaseqv2_20130521_brca_rank_sample_glm.csv")

glmminimum<-apply(glmres, 1, function(x){
  return (min(x[2:6]))
})

adjustpvalue<-p.adjust(glmminimum, method="bonferroni")
glmsig<-adjustpvalue < 0.01
glmpre<-cbind(glmres, glmminimum, adjustpvalue)
glmpost<-glmpre[glmsig,]
glmrank<-order(glmpost[,7])
glmfinal<-glmpost[glmrank,]

