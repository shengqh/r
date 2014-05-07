setwd("H:/shengquanhu/projects/tuchengjian/20140217_labelfree/")

data<-read.csv("map_linear_regression.csv", header=T)

proteins<-data$Protein.accession.numbers

mapdata<-data.frame(MAP_a=apply(data,1,function(x){  sum(as.numeric(x[2:6])) }),
                    MAP_b=apply(data,1,function(x){  sum(as.numeric(x[7:11])) }))

uniqpro<-unique(proteins)
ret<-data.frame(Protein=uniqpro, 
                Ratio=rep(1, length(uniqpro)), 
                StdErr=rep(0, length(uniqpro)), 
                tValue=rep(0, length(uniqpro)), 
                pValue=rep(1, length(uniqpro)), 
                PepCount=rep(0, length(uniqpro)))

index<-3
for(index in c(1:length(uniqpro))){
  pro<-uniqpro[index]
  prodata<-mapdata[proteins == pro,]
  rl<-lm(MAP_b~-1+MAP_a, data=prodata)
  srl<-summary(rl)
  coeff<-srl$coefficients
  ret[index,2:6]<-c(coeff[1], coeff[2], coeff[3], coeff[4], nrow(prodata))
}

write.csv(ret, "map_linear_regression_result.csv", row.names=F);
