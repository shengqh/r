#data<-read.table("H:/shengquanhu/projects/somaticmutation/all_rsmc/result/TCGA-A7-A0D9-RNA-TP-NT/TCGA-A7-A0D9-RNA-TP-NT.bases",sep="\t",header=T)
#file<-"1_120284424_C_C_A"

data<-read.table("H:/shengquanhu/projects/somaticmutation/all_rsmc/result/TCGA-BH-A0C0-DNA-TP-NB/TCGA-BH-A0C0-DNA-TP-NB.bases",sep="\t",header=T)
file<-"10_43597801_C_C_A"

library("brglm")
library("stringr")


parts<-unlist(strsplit(file,'_'))
chr<-parts[1]
locus<-parts[2]
ref<-parts[3]
major<-parts[4]
minor<-parts[5]

filedata<-data[data$Identity == file,]
filedata$SAMPLE<-factor(filedata$SAMPLE,levels=c("NORMAL", "TUMOR"))
filedata$Strand<-as.factor(filedata$Strand)
filedata$Base<-factor(filedata$Base,levels=c(major,minor))
filedata$Position<-as.factor(filedata$Position)

majordata<-filedata[filedata$Base == major,]
minordata<-filedata[filedata$Base == minor,]
hasScore = mean(minordata$Score) < mean(majordata$Score)

tumordata<-filedata[filedata$SAMPLE=="TUMOR",]
hasStrand = length(unique(tumordata$Strand)) > 1
hasPosition = length(unique(tumordata$Position)) > 1

filter<-"PASS"

tb<-table(filedata$Base, filedata$SAMPLE)

pvalue.fisher.group<-fisher.test(tb)$p.value

pvalue.brglm.group<-""
pvalue.brglm.strand<-""
pvalue.brglm.position<-""
pvalue.brglm.score<-""

baseformula<-"Base ~ SAMPLE";
res <- try(fit0<-brglm(baseformula, family=binomial, data=filedata))
addfake<-FALSE
failed<-FALSE
if(class(res) == "try-error" || !fit0$converged){
  if(tb[minor, "NORMAL"] == 0){
    normaldata<-filedata[filedata$SAMPLE=="NORMAL",]
    
    fake<-minordata[1,]
    fake$SAMPLE<-"NORMAL"
    fake$Score<-round(mean(minordata$Score))
    tbstrand<-table(normaldata$Strand)
    fake$Strand<-ifelse(tbstrand["FORWARD"] < tbstrand["REVERSE"], "FORWARD", "REVERSE")
    fake$Position<-"MIDDLE"      
    filedata<-rbind(filedata, fake)
    print ("Add fake normal data.")
    addfake<-TRUE
  }
  
  if(tb[major, "TUMOR"] == 0){
    fake<-majordata[1,]
    fake$SAMPLE<-"TUMOR"
    fake$Score<-round(mean(majordata$Score))
    tbstrand<-table(tumordata$Strand)
    fake$Strand<-ifelse(tbstrand["FORWARD"] < tbstrand["REVERSE"], "FORWARD", "REVERSE")
    fake$Position<-"MIDDLE"      
    filedata<-rbind(filedata, fake)
    print ("Add fake tumor data.")
    addfake<-TRUE
  }
  
  if(addfake){
    res <- try(fit0<-brglm(baseformula, family=binomial, data=filedata))
    if(class(res) == "try-error"){
      failed<-TRUE
    }
  }else{
    failed<-TRUE
  }
}

if(!failed){
  if(fit0$converged && coef(summary(fit0))[2, 4] < pvalue){
    formulas<-c("Base ~ SAMPLE + Score + Strand + Position",
                "Base ~ SAMPLE + Strand + Position",
                "Base ~ SAMPLE + Score + Position",
                "Base ~ SAMPLE + Score + Strand",
                "Base ~ SAMPLE + Position",
                "Base ~ SAMPLE + Strand",
                "Base ~ SAMPLE + Score",
                "Base ~ SAMPLE")
    
    if(!hasScore){
      formulas<-formulas[!str_detect(formulas, "Score")]
    }
    
    if(!hasStrand){
      formulas<-formulas[!str_detect(formulas, "Strand")]
    }
    
    if(!hasPosition){
      formulas<-formulas[!str_detect(formulas, "Position")]
    }
    
    for(curformula in formulas){
      cat("\t", curformula, "\n")
      
      if(curformula == baseformula){
        fit<-fit0
        break
      }
      
      res <- try(fit<-brglm(curformula, family=binomial, data=filedata))
      if(class(res) != "try-error" && fit$converged){
        break
      }
    }
  }else{
    cat("\t", baseformula, "\n")
    fit<-fit0
  }
}

pvalue.brglm.converged<-ifelse(failed, FALSE, fit$converged)

if(pvalue.brglm.converged){
  fit.coef=coef(summary(fit))
  for(index in c(1:nrow(fit.coef))){
    if(rownames(fit.coef)[index] == "SAMPLETUMOR"){
      pvalue.brglm.group<-fit.coef[index, 4] 
      pvalue.brglm_logistiX_group<-pvalue.brglm.group
    }else if(rownames(fit.coef)[index] == "Score"){
      pvalue.brglm.score<-fit.coef[index, 4] 
    }else if(rownames(fit.coef)[index] == "StrandREVERSE"){
      pvalue.brglm.strand<-fit.coef[index, 4] 
    }else if(rownames(fit.coef)[index] == "PositionTERMINAL"){
      pvalue.brglm.position<-fit.coef[index, 4] 
    }
  }
}

filter<-ifelse(failed, "GLM_FAILED", ifelse(pvalue.brglm.converged, "PASS", "NOT_CONVERGED"))

v<-c(chr, 
     locus, 
     locus, 
     major, 
     minor, 
     ref,
     tb[major, "NORMAL"], 
     tb[minor, "NORMAL"],
     tb[major, "TUMOR"], 
     tb[minor, "TUMOR"],
     pvalue.fisher.group, 
     pvalue.brglm.converged,
     pvalue.brglm.group, 
     pvalue.brglm.score, 
     pvalue.brglm.strand, 
     pvalue.brglm.position, 
     NA,
     filter, file)
