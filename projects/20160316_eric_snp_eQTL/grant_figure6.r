setwd("Y:/shengq1/20160316_eric_snp_eQTL")
#setwd("H:/shengquanhu/projects/20160316_eric_snp_eQTL")

require(ggplot2)
require(qqman)
require(MASS)
require(rms)

snp_genotype_file_name = "preliminary_filtered.genotype.txt"
bim_file_name="Hamid_HCE_Project_PLINK_271013_0138/preliminary_filtered.bim"
gene_expression_file_name<-"preliminary_filtered.gene_expression.txt"

snp_data<-read.table(snp_genotype_file_name, sep="\t", header=T, check.names=F, row.names=1)
snpbim<-read.delim(bim_file_name, header=F, row.names=2)
gene_expression<-read.table(gene_expression_file_name, sep="\t", header=T, check.names = F, row.names=1)

sorted_samples=sort(colnames(snp_data))
snp_data=snp_data[,sorted_samples]
gene_expression=gene_expression[,sorted_samples]

genesymbol<-"PTGIR"

cursnps=data.frame(snp=c("exm2268201", "rs6509226", "rs99595"), fdr=c(0.0461, 0.0262, 0.0238), stringsAsFactors = F)
curgene=genesymbol
cursnp=cursnps$snp[1]
associationFDR=cursnps$fdr[1]
gene_expression<-log(gene_expression)

alldat=NA
gtnames=c()
for(idx in c(1:3)){
  cursnp<-as.character(cursnps[idx, "snp"])
  associationFDR<-as.numeric(cursnps[idx, "fdr"])
  
  alleles<-snpbim[cursnp,]
  
  cat(cursnp, "\t", curgene, "\n")
  
  gt<-as.character(snp_data[cursnp,])
  ge<-as.numeric(gene_expression[curgene,])
  gtna<-is.na(gt) | gt=="NA"
  gt<-gt[!gtna]
  ge<-ge[!gtna]
  
  g0<-paste0(alleles$V6, alleles$V6)
  g1<-paste0(alleles$V6, alleles$V5)
  g2<-paste0(alleles$V5, alleles$V5)
  
  g0count<-paste0(g0, "(N=",length(gt[gt=="0"]), ")")
  g1count<-paste0(g1, "(N=",length(gt[gt=="1"]), ")")
  g2count<-paste0(g2, "(N=",length(gt[gt=="2"]), ")")
  gtnames=c(gtnames, g0count, g1count, g2count)

  dat<-data.frame(GenoType=gt, GeneExpression=ge, stringsAsFactors = F)
  dat$SNP=cursnp
  
  if(is.na(alldat)){
    alldat=dat
  }else{
    alldat=rbind(alldat, dat)
  }
}

png(file=paste0("figure6b.png"), height=900, width=2700, res=300)
g<-ggplot(alldat, aes(GenoType, GeneExpression)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.5)) +
  facet_grid(.~SNP) +
  xlab("Genotypes") +
  ylab("PTGIR Expression") +
  theme_grey() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

snps=unique(alldat$SNP)
ann_text1 <- data.frame(GenoType = rep("1",6),
                        GeneExpression = c(rep(-1.5,3), rep(-6.5,3)),
                        lab=c("P=.0461","P=.0262","P=.0238",paste0("Genotype of ", snps) ),
                        SNP=rep(snps,2))
ann_text2 <- data.frame(GenoType = rep(c("0","1","2"),3),
                        GeneExpression = rep(-6,9),
                        lab=gtnames,
                        SNP=rep(snps,1, each=3))

ann_text = rbind(ann_text1, ann_text2)

g=g+geom_text(data = ann_text,aes(label =lab) )
print(g)
dev.off()
