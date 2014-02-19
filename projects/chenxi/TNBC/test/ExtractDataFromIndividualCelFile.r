library ("affy")
library("Biobase")

getCelData<-function(celfile){
  expfile<- paste0(celfile,".tsv",sep="");
  
  if (!file.exists(expfile)){
    cat("Reading ", celfile, " ...\n", sep="");
    
    #read cel file
    data<-ReadAffy(filenames=celfile);
    
    rmadata<-rma(data,normalize=FALSE,background=TRUE);
    
    cat("Writing ", expfile, " ...\n", sep="");
    
    write.exprs(rmadata,file=expfile);
    
    #clean data
    rm(rmadata);
    rm(data);
  }
}

affymetrix_get_data<-function(curdir){
  celfiles<-list.celfiles(curdir,full.names=1);
  
  if(length(celfiles) == 0){
    return(1);
  }
  
  cels<-data.frame(celfiles);
  
  apply(cels,1,getCelData);
}

rootdir<-'I:/projects/BreastCancer/Dataset';

dirs<-data.frame(list.dirs(path=rootdir));

apply(dirs,1, affymetrix_get_data);


