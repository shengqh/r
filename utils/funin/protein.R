#normalization

x<-read.delim('clipboard',header=FALSE)
x=as.matrix(x)
x_norm<-normalize.quantiles(x)
write.table(x_norm,file="D:/x_norm.xls",quote=F,col.names=F,row.names=F,append=T)
save(x_norm,file="D:/raw.log2.norm.RData")
pdf("D:/raw.log2.norm.pdf")
boxplot(data.frame(x_norm))
dev.off()

#anova_adjst(BH)

norm<-read.delim('clipboard',header=FALSE)
norm=as.matrix(norm)

 

 tissue<-c(rep("CB",3),rep("PFC",3),rep("LV",3))
 aov_function<-function(x,t)
 {
 datas<-data.frame(gene=x,tissue=factor(t))
 d<-aov(gene~tissue,datas)
 result<-summary(d)[[1]][1,5]
 }
 
aov_results<-apply(norm,1,aov_function,t=tissue)
results=as.matrix(aov_results)
q=p.adjust(results,method="BH")
write.table(results,file="D:/aov_results.tsv",quote=F,col.names=F,row.names=F,append=T)
write.table(q,file="D:/aov_results_adj.tsv",quote=F,col.names=F,row.names=F,append=T)
diff_gene<-names(q[q<=0.01])
write.table(diff_gene,file="D:/different_sig_name.tsv",quote=F,col.names=F,row.names=F,append=T)


#PCA_method 2
x<-read.delim('clipboard',header=FALSE)
fix(x)
y<-as.matrix(x)
fix(y)
plot(prcomp(y))
y1<-prcomp(t(y))
plot(y1$x[,1],y1$x[,2],col=c('red','red','red','blue','blue','blue','green','green','green'),cex=3)
text(y1$x[,1]+2,y1$x[,2]+2,c("CB1","CB2","CB3","PFC1","PFC2","PFC3","LV1","LV2","LV3"),cex=0.8)

#merge files
 x<-read.delim('clipboard',header=TRUE)
> x=data.frame(x)
> fix(x)
>  y<-read.delim('clipboard',header=TRUE)
> y=data.frame(y)
> fix(y)
> z=merge(x,y,all=TRUE,sort=FALSE)
> fix(z)
> write.table(z,file="D:/merged.tsv",quote=F,col.names=F,row.names=F,append=T)

#p-value_adjust
p=read.delim('clipboard',header=FALSE)
> p=as.matrix(p)
> p1=p.adjust(p,"BH")
> write.table(p1,file="D:/p1.tsv")

#fourfold table(fisher_test)

data=read.delim('clipboard',header=FALSE)
data=as.matrix(data)
sig=apply(data,2,function(x)sum(data[which(x==1),1]))
number=apply(data,2,sum)
sig=matrix(sig)
number=matrix(number)
fix(sig)
fix(number)
combine=cbind(sig,number)
colnames(combine)=c("diff_num","target_num")
write.table(combine,file="D:/merged.tsv",sep="\t")

#arrange or calculate manually based on the combine table


tb=read.delim('clipboard',header=FALSE)
tb=as.matrix(tb)
fix(tb)
for(i in 1:nrow(tb)){
b=matrix(c(tb[i,1],tb[i,3],tb[i,2],tb[i,4]),nrow=2,ncol=2)
p.value<-fisher.test(b)$p.value
write.table(p.value,file="D:/p.value.tsv",append=T,col.names=F,row.names=F)}