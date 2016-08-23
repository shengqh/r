setwd("H:/shengquanhu/projects/20160803_anna_proteomics_metabolomics")

source("E:/sqh/programs/r/projects/20160803_anna/metabolomics_data_preparation.r")

library(glmnet)

varselect<-function(datamat, outcome, alpha){
  
  
  set.seed(2016)
  
  fit<-glmnet(x=datamat, y=outcome, family="binomial", alpha=alpha) # lasso or elastic net
  
  nsteps <- 100
  b1<-coef(fit)[-1, 1:nsteps]
  matplot(1:nsteps, t(b1), type = "o", pch = ".", col = 4, xlab = "Step", ylab = "Coefficients", lty = 1)
  abline(h = 0, lty = 2, col=2)
  
  
  
  # cross-validation
  # sample size large enough; 10-fold CV
  fit.cv<-cv.glmnet(x=datamat, y=outcome, family="binomial", type="class", nfolds=10)       
  plot(fit.cv)
  title("10-fold cross-validation", line=2.5)
  
  #cvm:   The mean cross-validated error - a vector of length length(lambda).
  #lambda.min:   value of lambda that gives minimum cvm.
  #lambda.1se: largest value of lambda such that error is within 1 standard error of the minimum.
  
  fit.cv$cvm
  min(fit.cv$cvm)
  fit.cv$lambda.min    
  fit.cv$lambda.1se    
  
  # position of smallest cvm:
  winner<-order(fit.cv$cvm)[1]
  
  betas<-fit$beta[, winner]
  winners<-betas[betas!=0]
  
  print(winners)
  
  write.csv(winners, paste(deparse(substitute(outcome)), "_winners.csv"))
  
}


# run lasso
varselect(data.matrix(d4[,16:184]), d4$PVH_PAH_lump, alpha=1)
varselect(data.matrix(d5[,16:184]), d5$PAH_NoPH, alpha=1)
varselect(data.matrix(d6[,16:184]), d6$PVH_Combined, alpha=1)
varselect(data.matrix(d7[,16:184]), d7$PAH_Combined, alpha=1)
varselect(data.matrix(d8[,16:184]), d8$PAH_PVH, alpha=1)
