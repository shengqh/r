setwd("I:\\projects\\GRO-Seq\\data\\range");

p<-read.table("mtg8_chr3121886990122714378.p.basecount",header=F,sep="\t");
r<-read.table("mtg8_chr3121886990122714378.r.basecount",header=F,sep="\t");

png(filename="mtg8_chr3121886990122714378.png",width=4000,height=3000,units="px",res=300);
plot(p[,2],p[,4],type="h",col="red",xlab="Base",ylab="Read Count",ylim=c(-50,50));
lines(r[,2],-r[,4],type="h", col="blue");
dev.off();




