###########################


load("/Users/YingHuang/data/415 final/FinalQ1.Rdata")
head(FinalQ1)

rm(list=ls())


##################
# 1.  One-liners:  Write R code to do the following tasks.  For full credit, 
# each must be coded in a single line of R code, fewer than 80 characters long 
# (code such as function(x) immediately followed by text is okay).  Of course, 
# simply printing out the results is right out.  The data set FinalQ1.Rdata, 
# available on the Blackboard site, contains example data which you can use 
# to test your code.  Please do not submit the output of your R code-
# just the code itself.  (30 pts total, 3 each)

# j.  Given vector yval of values of an outcome variable, and a matrix X 
# containing 5 columns of predictor variables, report the coefficients from
# an ordinary linear regression of y on each possible pair of columns of x.  
# That is, find the coefficients from linear models (ten in all) of  
# yval ~ X[,1] + X[,2],  yval ~ X[,1] + X[,3], ..., yval ~ X[,4] + X[,5].
unlist(lapply(c(1:4), function(i){unlist(lapply(c((i+1):5), function(j){lm(yval~x[,i]+x[,j])$}))}))

# k.  Bonus (3 extra points):  Given a vector xint of integers between
# 1 and 9999, e.g., xint = c(27, 181, 1009, 46, 2004, 7), generate a vector y 
# of four-digit numbers with padded zeros, e.g., 
# ypad = c("0027", "0181", "1009", "0046", "2004", "0007"), 
# in fewer than 30 characters of code and without using the nchar function.
ypad = formatC(xint, width = 4, format = "d", flag = "0")

########################
# 2.  Simulation study (30 pts) 
# 2.a and 2.b
rm(list=ls())
library(MASS)
in.ci <- function (x, ci) x > ci[1] && x < ci[2]

nrep <- 5000 # Number of replicates
a <- -2 # Alpha and beta
b <- 3
n <- 5 # Sample size; 

conf <- matrix(0,nrow=3,ncol=2)  # Matrix for 95% CI coverage probability for each error type
for (i in 1:nrep) {
  errmat <- matrix(c(mvrnorm(n,rep(0,5), 4*diag(5)),
                     mvrnorm(n,rep(0,5), 4*matrix(0.8,5,5)+diag(0.2,5)),
                     mvrnorm(n,rep(0,5), 4*0.8^(outer(1:5,1:5,function(x,y) abs(x-y))))),
                   nrow=n)

  x <- rnorm(n,mean=10,sd=3)
  y <- matrix(a + b*x + errmat, nrow=n)
  for(j in 1:3) {
    ci <- confint(lm(y[,j]~x))
    if(in.ci(a,ci[1,])) conf[j,1] <- conf[j,1] + 1
    if(in.ci(b,ci[2,])) conf[j,2] <- conf[j,2] + 1
  }
  if(round(i/100) == i/100) cat(i,"\n")
}

# Create table of 95% coverage probabilities
conf
conf <- conf/nrep
dimnames(conf) <- list(c("identity","exchangeable","AR1"), c("alpha","beta"))
conf

# output
# alpha   beta
# identity     0.9532 0.9548
# exchangeable 0.9512 0.9502
# AR1          0.9490 0.9494

# 2.c




########################
# 3.  Normal approximations (20 pts)

try(quartz())
par(mfrow=c(3,3), mai=rep(0.07,4), omi=c(1,1,0,0), las=1, 
    cex.lab=1.4, cex.axis=1.4, ljoin=1, lend=1)
for(i in 1:3){
  lam<-c(5,15,30)
  x <- seq(0:100)
  hx <- dnorm(x,mean=lam[i],sd=(lam[i])^(1/2))
  plot(x, hx, type="l", lty=1,lwd=3, col="red",xlab="x value",
       ylab="Density", main=paste0("Pois","(",lam[i],")"))
  lines(x, dpois(x,lam[i]), lwd=3, col="blue") 
}

for(i in 4:6){
  n<-c(25,50,100)
  p<- 0.8
  x<-seq(0,100,by=1)
  hx<-dnorm(x,mean=n[i]*p,sd=(n[i]*p*(1-p))^(1/2))
  plot(x, hx, type="l", lty=1,lwd=3, col="red",xlab="x value",
       ylab="Density", main=paste0("Binom","(",n[i],")"))
  lines(x, dbinom(x,n[i],p), lwd=3, col="blue") 
}

for(i in 7:9){
  df<-c(30,100,300)
  x <- seq(0:100)
  hx <- dnorm(x,mean=df[i],sd=(2*df[i])^(1/2))
  plot(x, hx, type="l", lty=1,lwd=3, col="red",xlab="x value",
       ylab="Density", main=paste0("Chisq","(",df[i],")"))
  lines(x, dchisq(x,df[i]), lwd=3, col="blue") 
}


