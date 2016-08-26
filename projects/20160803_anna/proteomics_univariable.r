setwd("H:/shengquanhu/projects/20160803_anna_proteomics_metabolomics")

source("E:/sqh/programs/r/projects/20160803_anna/proteomics_data_preparation.r")

#====================================================================
## univariable

library(rms)

## multiple imputation for age and BMI

a<-aregImpute(~Age + Gender + BMI + NoPH + PVH + Combined + PAH, n.impute=20, nk=0, data=d2) 

####################################################

fit.MI<-function(data, outcome, var){
  
  f<-fit.mult.impute(outcome ~ Age + Gender + BMI + var, fitter=lrm, xtrans=a, data=data, pr=F)
  
  #f<-lrm(outcome ~ d2$Age + d2$Gender + d2$BMI + var)
  #f<-glm(outcome ~ d2$Age + d2$Gender + d2$BMI + var, family=binomial(link='logit'))
  #print(f)
  return(f)
}


fit<-function(data, outcome, var){
  
  f<-summary(glm(outcome ~ data$Age + data$Gender + data$BMI + var, family=binomial(link='logit')))
  #print(f)
  return(f)
}

#==============================================================================

comp<-function(data, outcome, filename){
  coefficient<-NA
  pvalue<-NA
  
  for (i in 11:ncol(data)){
    
    output<-try(fit.MI(data, outcome, data[,i]))
    #print(output)
    
    if (is(output, "try-error")) {
      coefficient[i-10]<-NA
      pvalue[i-10]<-NA
    } 
    else 
    {  
      coefficient[i-10]<-output$coef[5]
      pvalue[i-10]<-anova(output)[4,3]   # multiple imputation
    } 
  }
  
  pvalue.adj<-p.adjust(pvalue, method="BH")
  res<-data.frame(proteome=colnames(data)[11:ncol(data)], coefficient, pvalue, pvalue.adj)
  res<-res[order(res$pvalue),]
  write.csv(res, filename, row.names=F)
  
  res<-res[!is.na(res$coefficient),]

  rnk<-data.frame(Protein=res$proteome, SignedPvalue=-log(res$pvalue) * sign(res$coefficient))
  write.table(rnk, paste0(filename, ".rnk"), col.names = F, row.names=F, sep="\t", quote=F)
  
}

comp(d4, d4$PVH_PAH, "proteomics_univariable_PVH_vs_PAH_with_imputation.csv")
comp(d5, d5$PAH_NoPH, "proteomics_univariable_PAH_vs_NoPH_with_imputation.csv")
comp(d6, d6$PVH_Combined, "proteomics_univariable_PVH_vs_Combined_with_imputation.csv")
comp(d7, d7$PAH_Combined, "proteomics_univariable_PAH_vs_Combined_with_imputation.csv")

# reference groups:

contrasts(as.factor(d4$PVH_PAH))
contrasts(as.factor(d5$PAH_NoPH))
contrasts(as.factor(d6$PVH_Combined))
contrasts(as.factor(d7$PAH_Combined))
