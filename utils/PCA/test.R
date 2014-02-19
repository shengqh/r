#引用函数文件
source("PCA_zhaoshilin.R");

#读取文件，按Tab分隔，第一列为row name，第一行为col name
data <- read.table("IMDL_data_662protein.txt",header=T,sep="\t",row.names=1)

#调用函数，前17个为一个状态，其余为第二个状态
PCA_2(data, number=17);
