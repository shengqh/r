root<-"d:/tmp";

fpkmfiles<-list.files(path=root, pattern="genes.fpkm_tracking", all.files=TRUE, full.names=TRUE,recursive=TRUE)

mybreaks = c(-100,0.5 * (-10:40),100);

maxcount<-0;

legends<-c();

colors<-topo.colors(length(fpkmfiles));

h<-list();

n=0;
for (file in fpkmfiles){
  n=n+1;

  genes<-read.table(file, sep="\t", header=1, row.names=1);

  fpkm<-as.matrix(genes["FPKM"]);

  h[[n]]<-hist(log(fpkm),breaks=mybreaks,plot=FALSE);
  
  maxcount<-max(maxcount, h[[n]]$counts);

  legends<-c(legends, basename(dirname(file)));
}

png(filename=paste0(root, "/fpkm.png"), width=4000, height=3000, res=300);

plot(0,xlim=c(-5,20),ylim=c(0,maxcount),type="n",xlab="log(FPKM)",ylab="Count",main="FPKM QC");

n=0;
for (file in fpkmfiles){
  n = n+1;
  
  hn<-h[[n]];
  
  counts<-hn$counts;
  
  y<-max(counts);
  
  xx<-hn$breaks[2:length(mybreaks)];
  
  lines(xx,hn$counts,col=colors[n]);
  
  kk<-c(1:(length(mybreaks)-1));
  
  for(k in kk){
    if (counts[k] == y){
      text(xx[k], y, legends[n]);
      break;
    }
  }
}

dev.off();
