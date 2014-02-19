# TODO: Add comment
# 
# Author: sqh
###############################################################################


inputfilename = "E:/sqh/Science/Project/Shift/10ppm/modif.txt";
g<-read.table(inputfilename ,header=T,sep="\t");
charge1<-subset(g,Charge==1);
show(nrow(charge1));
charge2<-subset(g,Charge==2);
show(nrow(charge2));
charge3<-subset(g,Charge>=3);
show(nrow(charge3));

