mycluster=function(filename, reverse=F, logdata=T, scale=F, jpeg=F, maintitle="", xlabel="", clusterCount=2){
data <- read.table(filename,header=T,sep="\t",row.names=1)

if(reverse==T){	data<-t(data); }

logname="";
if(logdata==T){ 
	logname = ".log";
	Data<-log(data); 
}

if(scale==T){data<-scale(data);}

# ��ɫ�ռ�
library(colorspace)
# ��ʾ���ݼ��Ľṹ
str(data)

# ϵͳ����
# �����һЩ��Ҫ�ĺ���
library(cluster)
library(rattle)
#ϵͳ���ຯ���ڰ�amap��
require(amap, quietly=TRUE)
#�������а�fpc�ṩ
require(fpc, quietly=TRUE)
#��ͼ ��cba��
require(cba, quietly=TRUE)

chcluster <- hclusterpar(na.omit(data), method="manhattan", link="ward", nbproc=2)

# ��������
centers.hclust(na.omit(data), chcluster, clusterCount)

#��������ͼ �þ�����ʾ������

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
