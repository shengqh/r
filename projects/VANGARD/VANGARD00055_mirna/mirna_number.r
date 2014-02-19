dir<-"H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna_v2/summary_1mm/result/"
files<-list.files(dir, pattern=".info$", full.name=T)
file<-files[1]

finaltable<-NA
for(file in files){
  data<-read.table(file, sep="\t", header=T, row.names=1, check.names=F)
  data$Dataset<-file_path_sans_ext(basename(file))
  
  if(is.na(finaltable)){
    finaltable<-data
  }else{
    finaltable<-rbind(finaltable, data)
  }
}

write.csv(finaltable, "H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna_v2/summary_1mm/summary.csv")
