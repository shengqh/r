#part1:functions

# 1. Read the example and modify it so
# that instead of finding the root of x^3 + 2 * x^2 - 7 you find the
# minimum of f = (x-3)^4 + 7*(x-2)^2, i.e. the point where df/dx = 0.  

#Newton’s method for root finding Example 4.7
# Suppose f(x) = x3 +2x2 −7. Then, if x0 is close enough to one of the three roots of this equation,
# xn = xn−1 − xn3−1 +2xn2−1 −7 3x2 + 4x
# will converge to a root.
# An R version of this could be implemented as follows:
# > x <- x0
# > f <- xˆ3 + 2 * xˆ2 - 7
# > tolerance <- 0.000001
# > while (abs(f) > tolerance) {
#   + f.prime<-3*xˆ2+4*x + x<-x-f/f.prime
#   + f<-xˆ3+2*xˆ2-7
#   +}
# >x

###let f = df/dx = 4 * (x-3)^3 + 14*(x-2)
###let initial x is 4
x0<-4
x <- x0
f <- 4 * (x-3)^3 + 14*(x-2)
tolerance <- 0.000001
while (abs(f) > tolerance) {
  f.prime<-9*(x-2)^2 +14
  x<-x-f/f.prime
  f <- 4 * (x-3)^3 + 14*(x-2)
}
x

# 2. Write a function mkPeople such that mkPeople(10) produces the list people 
# that problem 8(please see below) generated, but with a class attribute set to "people". 
# Write the function such that if no argument is present, a list of length 1 is
# produced. Use mkPeople() to make a list called folks of two people. 

# problem 8
#Create a list people of length 10 whose elements are named "person1" to 
# "person10"  where each element of the list is a list with 2 components: 
# id and age, filled, respectively, with IDs from above and
# # randomly sampled integers between 40 and 80.
# > people<-list()
# > for(i in c(1:10)){
#   > lst<-list(id=paste0("person", i), age=sample(c(40:80), 1))
#   > people[[i]]<-lst
#   > }

mkPeople<-function(number=1){
  people<-list()
  for(i in c(1:number)){
    name<-paste0("person", i)
    lst<-list(id=name, age=sample(c(40:80), 1))
    people[[name]]<-lst
  }
  return (structure(people, class="people"))
}
people1<-mkPeople()
people1
people2<-mkPeople(2)
people2

# 3. Write a function addBW such that addBW(folks) returns the input list
# with a new field called BW that holds rnorm(5,mean=rnorm(1,mean=70,sd=10),sd=2) 
# as the person's last 5 body weight measurements in kg. Use it to update folks.

addBW<-function(folks){
  for(i in c(1:length(folks))){
    folks[[i]]$BW<-rnorm(5,mean=rnorm(1,mean=70,sd=10),sd=2)
  }
  return (folks)
}
people3<-mkPeople(3)
people3bw<-addBW(people3)
people3bw


# 4. Write a function addBWmean such that addBWmean(people) returns the input 
# list with a new field called BWmean that holds the person's mean BW. Use it to
# update the folks output of problem 3.

addBWmean<-function(folks){
  for(i in c(1:length(folks))){
    folks[[i]]$BWmean<-mean(folks[[i]]$BW)
  }
  return (folks)
}
people3bwmean<-addBWmean(people3bw)
people3bwmean

# 5. Write an S3 summary method summary.people() to summarize the object. 
# The method should cat() the average population age and weight to the screen 
# and return a data.frame where the row names are "person1","person2", etc, and 
# the columns are the patient's ID, age, and average body weight rounded to
# one decimal point using round().

# The output of summary(folks) using the folks output of problem 4 should thus
# be similar to this:

# > summary(folks)
# The mean age of the population is 68.5 years.
# The mean body weight of the population is 84.9  kg.
# 
#            ID age meanBW
# person1 27AD3  74   91.6
# person2 AEBD5  63   78.2

summary.people<-function(folks){
  folklen<-length(folks)
  agesum<-0
  weightsum<-0
  for(i in c(1:folklen)){
    agesum<-agesum + folks[[i]]$age
    weightsum<-weightsum +folks[[i]]$BWmean
  }
  cat("The mean age of the population is", formatC(agesum / folklen, digits=1, format="f"), "years\n")
  cat("The mean body weight of the population is", formatC(weightsum / folklen,digits=1, format="f"), "years\n")
  cat("\n")
  cat("\tID\tage\tmeanBW\n")
  folknames<-names(folks)
  for(i in c(1:folklen)){
    cat(folknames[i], folks[[i]]$age, formatC(folks[[i]]$BWmean,digits=1, format="f") , "\n", sep='\t')
  }
}

summary(people3bwmean)

## part2:graphics. 
# For this problem set, please bundle any figures requested into a Word file and
# save it as a pdf. Both this pdf and R script are required. 

# 1. rgl is a 3D graphics package that allows you to rotate data to find
# best angles interactively. For example,
library(rgl)
bg3d("white")
with(iris,plot3d(Sepal.Width, Petal.Length,Petal.Width,
                 col=as.numeric(Species)))
#shows that three plant species in iris separate according to sepal and petal
#measurements

# study the rgl help pages to learn how to plot spheres instead of points. Find 
# a sphere radius size that you like and capture a good angle image using, for
# example, the Windows Accessories program Snipping Tool.

with(iris,plot3d(Sepal.Width, Petal.Length,Petal.Width, col=as.numeric(Species), type="s"))




# 2. One reason that ggplot2 cannot completely replace base graphics is that 
# many packages already have plot methods for specific classes of objects.
# Consider the following melanoma survival data and code. 

library(ISwR)
library(survival)
?melanom  # status = 1 = dead from cancer, 2= alive, 3= dead from other cause
head(melanom) # note that rows 5 and 6 are cases dead due to melanoma
(S=with(melanom,Surv(days,status==1))) 
# 10+ means we only know it would have been more than 10 days
dput(S) # object is class Surv, a matrix where 1 in second col=> dead with disease (dwd)
(kms=survfit(S~sex,data=melanom))
dput(kms) #class survfit
plot(kms) # plot method exists in survival package for survfit class objects 

# modify this plot to color code females (who have better apparent survival) as red and 
# males (who have worse apparent survival) as blue. Include a legend, axis labels, and a
# title. The legend should include the number of males and the number of
# females. Also, remove any excess white space that may exist on some of the
# margins.

par(mar=c(5,5,1,1))
plot(kms, col=c("red", "blue"), xlab="time", ylab="survival")
legend(100, .4, c("Female", "Male"), col=c("red", "blue"), lty=1) 

# 3. Study the code below 
load("/Users/radivot/CML.RData") # brings in dataframe d
head(d)
graphics.off()
(p <- ggplot(d,aes(x=age,y=incid,shape=Decade))+geom_point(size=4)
 + labs(title="CML Incidence",x="Age (years)",
        y=expression(paste("Cases per ",10^6," Person-Years")))    
 + scale_y_log10(limits=c(3,200)) )
(p=p + facet_grid(code ~ sex))
(p=p+theme(legend.position = c(.60, .86)) )

# modify this code to render the 3 decades in the colors red, green and blue