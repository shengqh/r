setwd("D:\\chungI\\Vandy\\Sequencing\\Somatic Mutation Data\\P2277_wsm")
fl<-list.files(, pattern=".wsm")
# Bias-reduced inferen
# install.packages("brglm")
# table(raw.data$Base)
# table( raw.data$Base, raw.data$SAMPLE)
library("brglm")
out<-matrix(NA, length(fl), 3)
#for (i in 1:3) {
for (i in 1:length(fl)) {
  raw.data<-read.table(fl[i], header=T)
  flag.start<-all( raw.data$IsStart == TRUE) | all( raw.data$IsStart == FALSE)
  flag.end<-all( raw.data$IsEnd == TRUE) | all( raw.data$IsEnd == FALSE)
  major<-names( which.max(table(raw.data$Base)) )
  minor<-names( which.min(table(raw.data$Base)) )
  raw.data$Base<-factor(raw.data$Base, levels=c(minor, major) )
  pvalue.fisher<-fisher.test( table( raw.data$Base, raw.data$SAMPLE))$p.value
  
  #fit<-glm(Base ~ factor(SAMPLE) + Score , family=binomial, data=raw.data)
  
  fit1<-brglm(Base ~ factor(SAMPLE) + Score , family=binomial, data=raw.data)
  pvalue.logistic.score<-coef(summary(fit1))[2, 4]
  
  if (flag.start | flag.end) {
    pvalue.logistic.all<-NA
  }else{
    fit2<-brglm(Base ~ factor(SAMPLE) + Score + factor(IsStart) + factor(IsEnd), family=binomial, data=raw.data)
    pvalue.logistic.all<-coef(summary(fit2))[2, 4]
  }
  out[i, ]<-cbind(pvalue.fisher, pvalue.logistic.score, pvalue.logistic.all)
}
colnames(out)<-c("Fisher", "logistic_score", "logistic_all")
rownames(out)<-fl
write.table(out, "Somatic_Mutation_Data")