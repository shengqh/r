library ("affy")
library ("affxparser")
library ("limma")
library ("hash")

affymatrix_get_chiptype<-function(curdir){
  celfiles<-list.celfiles(curdir,full.names=1);
  
  if(length(celfiles) == 0){
    return ('');
  }
  
  cels<-data.frame(celfiles);
  
  getChipType<-function(celfile){
    #read cel file
    data<-ReadAffy(filenames=celfile);
    
    #return cdfName
    chiptype<-cdfName(data);
    
    #clean data
    rm(data);
    
    #return chiptype
    return(chiptype);
  }
  
  return (getChipType(celfiles[1]));
}

rootdir<-'I:/projects/BreastCancer/Dataset';
dirs<-list.dirs(path=rootdir);

out<-t(sapply(1:2, function(i) c(dirs[i], affymatrix_get_chiptype(dirs[i])) ));
colnames(out)<-c('Dir','ChipType');

resfile<-paste0(rootdir, "/summary.csv");
write.csv(out, file=resfile, row.names=FALSE )



