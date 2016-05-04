setwd("Y:/shengq1/20160316_eric_snp_eQTL")

library(ggplot2)
library(qqman)
require(MASS)
require(rms)
require(biomaRt)
require(snpStats)

# adjust the limitation of y axis based on cis and trans pvalue
manhattanYlim<-15

# Distance for local gene-SNP pairs
cisDist = 1e6;

snp_genotype_file_names = c("preliminary_filtered.genotype.txt", "comprehensive_filtered.genotype.txt")
snps_location_file_names = c("preliminary_filtered.genotype.location.txt", "comprehensive_filtered.genotype.location.txt")
covariates_file_names=c("preliminary_filtered.covariates.all.txt","comprehensive_filtered.covariates.all.txt")
bim_file_names=c("Hamid_HCE_Project_PLINK_271013_0138/preliminary_filtered.bim","Hamid_HCE_Project_PLINK_271013_0138/comprehensive_filtered.bim")
gene_expression_file_names<-c("preliminary_filtered.gene_expression.txt", "comprehensive_filtered.gene_expression.txt")

drawpicture<-function(snpbim, filtered, name, snp_data, gene_data){
  dir.create(name, showWarnings = FALSE)
  
  cursnp<-as.character(filtered[1,1])
  curgene<-as.character(filtered[1,3])
  
  apply(filtered, 1, function(x){
    cursnp<-as.character(x[1])
    curgene<-as.character(x[3])
    associationFDR<-as.numeric(x[5])
    kruskalFDR<-as.numeric(x[8])
    
    alleles<-snpbim[cursnp,]
    
    cat(cursnp, "\t", curgene, "\n")
    
    gt<-as.character(snp_data[cursnp,])
    ge<-as.numeric(gene_data[curgene,])
    gtna<-is.na(gt)
    gt<-gt[!gtna]
    ge<-ge[!gtna]
    
    g0<-paste0(alleles$V5, alleles$V5)
    g1<-paste0(alleles$V5, alleles$V6)
    g2<-paste0(alleles$V6, alleles$V6)
    
    levels<-c()
    
    if(length(gt[gt=="0"]) > 0){
      g0count<-paste0(g0, "(N=",length(gt[gt=="0"]), ")")
      gt[gt=="0"]<-g0count
      levels<-c(levels, g0count)
    }
    
    if(length(gt[gt=="1"]) > 0){
      g1count<-paste0(g1, "(N=",length(gt[gt=="1"]), ")")
      gt[gt=="1"]<-g1count
      levels<-c(levels, g1count)
    }
    
    if(length(gt[gt=="2"]) > 0){
      g2count<-paste0(g2, "(N=",length(gt[gt=="2"]), ")")
      gt[gt=="2"]<-g2count
      levels<-c(levels, g2count)
    }
    
    dat<-data.frame(GenoType=gt, GeneExpression=ge)
    dat$GenoType<-factor(dat$GenoType, levels=levels)
    
    genename<-gsub("_", " ", curgene, perl=TRUE)
    
    png(file=paste0(name, "/", curgene, "_", cursnp, "_fdr", x[5], ".png"), height=2000, width=2000, res=300)
    g<-ggplot(dat, aes(GenoType, GeneExpression)) + 
      geom_boxplot(outlier.shape = NA) + 
      geom_point(position = position_jitter(width = 0.5)) +
      xlab(paste0("Genotype of ", cursnp)) + 
      ylab(paste0("Gene expression of ", genename)) +
      ggtitle(paste0("association.fdr=",format.pval(associationFDR), "; kruskal.fdr=",format.pval(kruskalFDR)))
    print(g)
    dev.off()
  })
}

do_eqtls<-function(snps, gene_covariants){
  result=snp_pos[,c(1:4)]
  result$pvalue=apply(snps, 1, function(x){
    gene_covariants$GenoType=x
    fit = lm( GenoType ~ GeneExpression + gender + family, data = gene_covariants );
    summary(fit)$coefficients["GeneExpression",4]
  })
  
  result$fdr=p.adjust(result$pvalue)
  result=result[order(result$pvalue),]
  return (result)
}

perform_eqtls<-function(eqtls_file_name_noext, snps, gene_covariants, fdr_threshold=0.05){
  allfile<-paste0(eqtls_file_name_noext, ".tsv")
  if(!file.exists(allfile)){
    cat('Calculating', eqtls_file_name_noext, 'eqtls\n')
    eqtls=do_eqtls(snps, gene_covariants)
    write.table(eqtls, file=allfile, sep="\t", row.names=F)
  }else{
    eqtls=read.delim(allfile, header=T)
  }
  
  kw_sig_file<-paste0(eqtls_file_name_noext, ".lm_sig.kw_sig.tsv")
  if(!file.exists(kw_sig_file)){
    lm_sig_file<-paste0(eqtls_file_name_noext, ".lm_sig.tsv")
    if(!file.exists(lm_sig_file)){
      sig_eqtls<-eqtls[eqtls$fdr <= fdr_threshold,]
      write.table(sig_eqtls, file=lm_sig_file, sep="\t", row.names=F)
    }else{
      sig_eqtls=read.delim(lm_sig_file, header=T)
    }
    
    sig_eqtls$kruskal = unlist(apply(sig_eqtls, 1, function(x){
      cursnp=x["snp"]
      gt<-snps[cursnp,]
      gtna<-is.na(gt)
      gt<-gt[!gtna]
      ge<-gene_covariants$GeneExpression[!gtna]
      kruskal.test(ge, as.factor(gt))$p.value
    }))    
    sig_eqtls$kruskalfdr=p.adjust(sig_eqtls$kruskal)
    kw_sig_eqtls<-sig_eqtls[sig_eqtls$kruskalfdr < fdr_threshold,]
    write.table(kw_sig_eqtls, file=kw_sig_file, sep="\t", row.names=F)
  }else{
    kw_sig_eqtls=read.delim(kw_sig_file, header=T)
  }
  
  list(eqtls=eqtls, kw_sig_eqtls=kw_sig_eqtls)
}

index=1
for(index in c(1:2)){
  bim_file_name<-bim_file_names[index]
  snpbim<-read.delim(bim_file_name, header=F, row.names=2)
  
  covariates_file_name = covariates_file_names[index]
  covariates_file_name_noext<-sub("[.][^.]*$", "", covariates_file_name, perl=TRUE)
  covariates<-read.delim(covariates_file_name, header=T, check.names = F, row.names=1)
  covariates<-t(covariates)
  
  gene_expression_file_name<-gene_expression_file_names[index]
  gene_expression_file_name_noext<-sub("[.][^.]*$", "", gene_expression_file_name, perl=TRUE)
  gene_expression<-read.table(gene_expression_file_name, sep="\t", header=T, check.names = F, row.names=1)
  genesymbols<-rownames(gene_expression)
  
  snp_genotype_file_name <-snp_genotype_file_names[index]
  snp_genotype_file_name_noext<-sub("[.][^.]*$", "", snp_genotype_file_name, perl=TRUE)
  snp_data<-read.table(snp_genotype_file_name, sep="\t", header=T, check.names=F, row.names=1)
  
  snps_location_file_name = snps_location_file_names[index]
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
  sorted_samples=sort(colnames(snp_data))
  
  snp_data=snp_data[,sorted_samples]
  gene_expression=gene_expression[,sorted_samples]
  covariates=covariates[sorted_samples,]
  
  gene_location_file_name<-paste0(gene_expression_file_name_noext, ".location.txt")
  if(!file.exists(gene_location_file_name)){
    gss<-sub("_.*$", "", genesymbols, perl=TRUE) 
    grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    genepos <- lapply(gss, function(x){
      getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
            filters="external_gene_name",
            values = x,
            mart=grch37)})
    genepos <- do.call("rbind", genepos)
    genepos$external_gene_name<-rownames(gene_expression)
    colnames(genepos)<-c("geneid",	"chr",	"s1",	"s2")
    write.table(genepos, file=gene_location_file_name, row.names = F, quote = F, sep="\t")
  }else{
    genepos = read.table(gene_location_file_name, sep="\t", header=T)
  }
  
  row.names(genepos)<-genepos$geneid
  genepos$cismin<-genepos$s1 - cisDist
  genepos$cismax<-genepos$s2 + cisDist
  
  genesymbol<-genesymbols[1]
  for(genesymbol in genesymbols){
    expression_file_name_noext<-paste0(gene_expression_file_name_noext, ".", genesymbol)
    cat(expression_file_name_noext, "\n")  
    
    genename<-gsub("_", " ", genesymbol, perl=TRUE)
    
    genedata<-unlist(gene_expression[genesymbol,,drop=T])
    geneloc<-genepos[genesymbol,]
    
    cursnppos<-cbind(data.frame(gene=genesymbol, geneloc=paste0(geneloc$chr, ":", geneloc$s1,"-",geneloc$s2), snploc=paste0(snpspos$chr, ":", snpspos$pos)), snpspos)
    cursnppos<-cursnppos[,c(4,3,1,2,5,6)]
    
    ge_covariants<-cbind(data.frame(GeneExpression=genedata), covariates)
    
    cis_snp_pos=cursnppos[which(cursnppos$chr==geneloc$chr[1]),]
    cis_snp_pos=cis_snp_pos[which(cis_snp_pos$pos>geneloc$cismin[1]),]
    cis_snp_pos=cis_snp_pos[which(cis_snp_pos$pos<geneloc$cismax[1]),]
    cis_snp = snp_data[cis_snp_pos$snp,]
    
    cis_result<-perform_eqtls(paste0(expression_file_name_noext, ".cis"), cis_snp, ge_covariants)
    cis_eqtls=cis_result$kw_sig_eqtls
    
    drawpicture(snpbim, cis_eqtls, "cis", snp_data, gene_expression)
    
    trans_snp_pos=cursnppos[! cursnppos$snp %in% cis_snp_pos$snp,]
    trans_result<-perform_eqtls(paste0(expression_file_name_noext, ".trans"), genesymbol, snp_data, trans_snp_pos, dat)
    trans_eqtls=trans_result$
    
    
    allsnps<-rbind(cis_eqtls, trans_eqtls)
    
    snppvalue<-data.frame(SNP=cis_eqtls$snp, CHR=snpchrs, BP=snpbps, P=allsnps)
    snppvalue<-snppvalue[snppvalue$CHR < 24,]
    snppvalue<-snppvalue[with(snppvalue, order(CHR, BP)),]
    
    png(file=paste0(expression_file_name_noext, ".", modelName, ".manhattan.png"), width=2000, height=2000, res=300)
    manhattan(snppvalue, ylim=c(0, manhattanYlim), main=genename)
    dev.off()
  }
}