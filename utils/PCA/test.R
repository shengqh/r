#���ú����ļ�
source("PCA_zhaoshilin.R");

#��ȡ�ļ�����Tab�ָ�����һ��Ϊrow name����һ��Ϊcol name
data <- read.table("IMDL_data_662protein.txt",header=T,sep="\t",row.names=1)

#���ú�����ǰ17��Ϊһ��״̬������Ϊ�ڶ���״̬
PCA_2(data, number=17);