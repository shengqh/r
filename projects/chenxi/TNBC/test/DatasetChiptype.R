library ("affy")
library ("affxparser")
library ("limma")


affymatrix_get_chip_type<-function(celfile){
  print(celfile);

  #result file
  resultfile<-paste0(celfile, ".tsv");
  
  if (file.exists(resultfile)){
    data<-readCelHeader(filename=celfile);
    
    chiptype<-data$chiptype;
    
    rm(data);
  }
  else{
    #read cel file
    data<-ReadAffy(filenames=celfile);
  
    #Normalizes data with 'rma' function and assigns them to ExpressionSet object
    eset<-rma(data);
  
    #save expression file
    write.exprs(eset, resultfile);
  
    #return cdfName
    chiptype<-cdfName(data);
  
    #clean data
    rm(data, eset);
  }
  
  #return chiptype
  return(chiptype);
}

affymatrix_count_chiptype<-function(curdir){
  print(curdir);
  
  celfiles<-data.frame(list.celfiles(curdir,full.names=1));
  
  if(nrow(celfiles) == 0){
    print(paste0('No cel file in directory ', curdir));
    return (data.frame());
  }

  chiptypes<-apply(celfiles, 1, affymatrix_get_chip_type);
  
  tc<-data.frame(table(chiptypes));
  colnames(tc)[1]<-'Chiptype';
  colnames(tc)[2]<-'Sample';
  
  dirname<-basename(curdir);
  dirnames<-data.frame(c(dirname, rep('',nrow(tc)-1)));
  colnames(dirnames)[1]<-'Dataset';
  
  tcc<-cbind(dirnames, tc);
  return(tcc);
}

rootdir<-'I:/projects/BreastCancer/Dataset';
#rootdir<-'I:/projects/BreastCancer/temp';

dirs<-list.dirs(path=rootdir,recursive=FALSE);

#cts<-apply(data.frame(dirs, 1, affymatrix_count_chiptype);

cts<-data.frame();
for(i in 1:length(dirs)){
  ct<-affymatrix_count_chiptype(dirs[i]);
  cts<-rbind(cts, ct);
}

write.csv(cts, file=paste0(rootdir,'/summary.csv'),row.names=FALSE);
