setwd("H:/shengquanhu/projects/20160803_anna_proteomics_metabolomics")

source("E:/sqh/programs/r/projects/20160803_anna/metabolomics_data_preparation.r")

fit.MI<-function(data, outcome, var){
  
  f<-fit.mult.impute(outcome ~ Age + Gender + BMI + var, fitter=lrm, xtrans=a, data=data)
  
  #f<-lrm(outcome ~ d2$Age + d2$Gender + d2$BMI + var)
  #f<-glm(outcome ~ d2$Age + d2$Gender + d2$BMI + var, family=binomial(link='logit'))
  #print(f)
  return(f)
}

comp<-function(data, outcome, filename){
  coefficient<-NA
  pvalue<-NA
  
  for (i in 16:184){
    
    output<-try(fit.MI(data, outcome, data[,i]))
    #print(output)
    
    if (is(output, "try-error")) {
      coefficient[i-15]<-NA
      pvalue[i-15]<-NA
    } 
    else 
    {  
      coefficient[i-15]<-output$coef[5]
      pvalue[i-15]<-anova(output)[4,3]   # multiple imputation
      #pvalue[i-15]<-output$coef[5,4]   # no imputation
    } 
  }
  
  pvalue.adj<-p.adjust(pvalue, method="BH")
  res<-data.frame(metabolome=colnames(data)[16:184], coefficient, pvalue, pvalue.adj)
  res<-res[order(res$pvalue),]
  print(res)

  write.csv(res, filename)
}

comp(d4, d4$PVH_PAH_lump, "PVH_vs_PAH_lump_univariable_results_with_imputation.csv")
comp(d5, d5$PAH_NoPH, "PAH_vs_NoPH_univariable_results_with_imputation.csv")
comp(d6, d6$PVH_Combined, "PVH_vs_Combined_univariable_results_with_imputation.csv")
comp(d7, d7$PAH_Combined, "PAH_vs_Combined_univariable_results_with_imputation.csv")
comp(d8, d8$PAH_PVH, "PAHnoCTD_vs_PVH_univariable_results_with_imputation.csv")

# reference groups:

contrasts(as.factor(d4$PVH_PAH_lump))
contrasts(as.factor(d5$PAH_NoPH))
contrasts(as.factor(d6$PVH_Combined))
contrasts(as.factor(d7$PAH_Combined))
contrasts(as.factor(d8$PAH_PVH))
