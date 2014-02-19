mycluster=function(filename, reverse=F, logdata=T, scale=F, jpeg=F, maintitle="", xlabel="", clusterCount=2){
data <- read.table(filename,header=T,sep="\t",row.names=1)

if(reverse==T){	data<-t(data); }

logname="";
if(logdata==T){ 
	logname = ".log";
	Data<-log(data); 
}

if(scale==T){data<-scale(data);}

# 彩色空间
library(colorspace)
# 显示数据集的结构
str(data)

# 系统聚类
# 聚类的一些必要的函数
library(cluster)
library(rattle)
#系统聚类函数在包amap中
require(amap, quietly=TRUE)
#聚类结果有包fpc提供
require(fpc, quietly=TRUE)
#绘图 需cba包
require(cba, quietly=TRUE)

chcluster <- hclusterpar(na.omit(data), method="manhattan", link="ward", nbproc=2)

# 聚类中心
centers.hclust(na.omit(data), chcluster, clusterCount)

#产生树形图 用矩形显示聚类结果

if(jpeg == T){
	outfilename = paste(filename, logname, ".cluster.jpg",sep="");
	jpeg(file=outfilename, width=480, height=640);
}

par(bg="grey")
plot(chcluster, main=maintitle, sub="", xlab=xlabel, hang=0)
rect.hclust(chcluster, k=clusterCount)

if(jpeg == T){
	dev.off();
}

}
