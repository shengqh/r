#this function is used to do fisher-exact-test for gene
groseq_fisherexact_pausing_forgene<-function(x, 
                                              tssCountIndex, 
                                              tssLengthIndex, 
                                              genebodyCountIndex, 
                                              genebodyLengthIndex){
  if (is.na(x[tssCountIndex]) || is.na(x[genebodyCountIndex]) || is.na(x[tssLengthIndex]) || is.na(x[genebodyLengthIndex])){
    result<-NA
  }
  else{
    c1 <- round(as.numeric(x[tssCountIndex]));
    c2 <- round(as.numeric(x[genebodyCountIndex]));
    l1 <- round(as.numeric(x[tssLengthIndex]));
    l2 <- round(as.numeric(x[genebodyLengthIndex]));
    
    if (c1 == 0 || c2 == 0 || l1 == 0 || l2 == 0){
      result<-NA
    }
    else{
      expectC1<-round(l1 * (c1 + c2) / (l1 + l2));
      expectC2<-(c1 + c2 - expectC1);
      v <- matrix(c(c1, expectC1, c2, expectC2), nrow=2);
      ft<-fisher.test(v);
      result<-ft$p.value;
    }
  }
}

#this function is used to do fisher-exact-test for dataset
groseq_fisherexact_pausing<-function(grodata, 
                                     tssCountIndex, 
                                     tssLengthIndex, 
                                     genebodyCountIndex, 
                                     genebodyLengthIndex, 
                                     newcolname){
  fs<-apply(grodata, 
                       1, 
                       groseq_fisherexact_pausing_forgene, 
                       tssCountIndex,
                       tssLengthIndex,
                       genebodyCountIndex,
                       genebodyLengthIndex);
  fs<-p.adjust(fs,"BH");
  ffs<-data.frame(fs);
  colnames(ffs) = newcolname;
  result<-cbind(grodata,ffs);
}

#this function is used to gene comparison between two conditions for gene
groseq_fisherexact_comparison_forgene<-function(x, 
                                                 tssCountIndex1, 
                                                 genebodyCountIndex1, 
                                                 tssCountIndex2, 
                                                 genebodyCountIndex2){
  if (is.na(x[tssCountIndex1]) || is.na(x[tssCountIndex2]) || is.na(x[genebodyCountIndex1]) || is.na(x[genebodyCountIndex2])){
    result<-NA
  }
  else{
    c1 <- round(as.numeric(x[tssCountIndex1]))
    g1 <- round(as.numeric(x[genebodyCountIndex1]))
    c2 <- round(as.numeric(x[tssCountIndex2]))
    g2 <- round(as.numeric(x[genebodyCountIndex2]))
    
    v <- matrix(c(c1, g1, c2, g2), nrow=2);
    ft<-fisher.test(v)
    result<-ft$p.value
  }
}

#this function is used to do fisher-exact-test for entire data
groseq_fisherexact_comparison<-function(grodata, 
                                        tssCountIndex1, 
                                        genebodyCountIndex1, 
                                        tssCountIndex2, 
                                        genebodyCountIndex2,
                                        newcolname){
  fs<-apply(grodata, 
                       1, 
                       groseq_fisherexact_comparison_forgene, 
                       tssCountIndex1, 
                       genebodyCountIndex1, 
                       tssCountIndex2, 
                       genebodyCountIndex2);
  fs<-p.adjust(fs,"BH");
  ffs<-data.frame(fs);
  colnames(ffs) = newcolname;
  result<-cbind(grodata,ffs);
}

library("stats");

#processing file
mydir<-"I:/projects/GRO-Seq/data/";

files<-c(
  'MTG8-gro-uniq-mapped_promoter50bp-count-forstat.csv',
  'groseq-3t3HDAC3-uniq.mapped-promoter50-count-forstat.csv'
);

isnorm<-c(
  0,
  0
);

for(i in 1:length(files)){
  currFile<-paste0(mydir,files[i]);
  isnormFile<-isnorm[i];
  
  mf<-read.table(currFile,header=T,sep=",")
  mf<-groseq_fisherexact_pausing(mf,2,3,6,7,'SH-1-Pausing-pValue');
  
  if(isnormFile){
    mf<-groseq_fisherexact_pausing(mf,12,13,17,18,'SH-2-Pausing-pValue');
    mf<-groseq_fisherexact_comparison(mf,2,6,12,17,'SH-1-2-Change-pValue');
  }
  else{
    mf<-groseq_fisherexact_pausing(mf,11,12,15,16,'SH-2-Pausing-pValue');
    mf<-groseq_fisherexact_comparison(mf,2,6,11,15,'SH-1-2-Change-pValue');
  }

  resultFile<-paste0(currFile,".pValue.csv");
  write.table(mf,file=resultFile,sep=",",row.names=F)
}