# TODO: Add comment
# 
# Author: sqh
###############################################################################



if(i == 1){
  screen(5)
  
  plot(0,xlim=c(0,4),ylim=c(0,1),type="n",xlab="log(Score)",ylab="Density")
  
  gFiltered<-subset(subset(g,Decoy=="True"), Modified=="True");
  lines(density(log(gFiltered[,1])), col="blue");
  
  gFiltered<-subset(subset(g,Decoy=="True"), Modified=="False");
  lines(density(log(gFiltered[,1])), col="brown");
  
  legend(1,1.0,legend= c(t2,t4),lty=4,col=c("blue","brown"));
}

if(i == 3){
  screen(6)
  
  plot(0,xlim=c(0,4),ylim=c(0,1),type="n",xlab="log(Score)",ylab="Density")
  
  gFiltered<-subset(subset(g,Decoy=="False"), Modified=="True");
  lines(density(log(gFiltered[,1])),col="red");
  
  gFiltered<-subset(subset(g,Decoy=="True"), Modified=="True");
  lines(density(log(gFiltered[,1])),col="blue");
  
  gFiltered<-subset(subset(g,Decoy=="False"), Modified=="False");
  lines(density(log(gFiltered[,1])),col="green");
  
  gFiltered<-subset(subset(g,Decoy=="True"), Modified=="False");
  lines(density(log(gFiltered[,1])),col="brown");
  
  legend(1,1.0,legend= c(t1,t2,t3,t4),lty=4,col=c("red","blue","green","brown"));
}
