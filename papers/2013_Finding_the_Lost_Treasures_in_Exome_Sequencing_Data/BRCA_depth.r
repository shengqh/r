data<-read.table("E:/sqh/Dropbox/papers/2013_Finding_the_Lost_Treasures_in_Exome_Sequencing_Data/BRCA_depth_matrix.txt", sep="\t", header=T, row.names=1)

mediandepth<-apply(data, 2, median)
ordered<-mediandepth[order(mediandepth)]
om<-as.matrix(ordered)
colnames(om)[1]<-"Median(Depth)"

write.csv(om, "E:/sqh/Dropbox/papers/2013_Finding_the_Lost_Treasures_in_Exome_Sequencing_Data/BRCA_median_depth.csv")

min(mediandepth)
max(mediandepth)