library(stats)
y<-c(-0.03141, 0.11393, 0.24809, 0.41024, 0.53651, 0.7031)
x<-c(0,8,24,48,72,120)
plot(y~x,type="l")
fit<-nls(y~a*(1-exp(-b*x)), start=list(a=0.8,b=0.01))
summary(fit)$parameters
