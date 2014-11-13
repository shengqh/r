
setwd("H:/shengquanhu/projects/chenxi/20140527_chenxi_pietenpol_p53")

chrs = c(1,3,6,7,9,10,11,12,15,16,17,18,19,20,22)
chrstrs=c("01","03","06","07","09","10","11","12","15","16","17","18","19","20","22")
pdf("STEP11_BeforeAndAfterPhasing.pdf")
for(index in c(1:length(chrs))){
  chr  =chrs[index]
  chrstr = chrstrs[index]

  prefile = paste0("STEP09_preimputation/Pietenpol_p53.", chr, ".gen.allele2frequency.tsv")
  postfile = paste0("STEP10_shapeit_gen/result/Pietenpol_p53.", chrstr, "/Pietenpol_p53.", chrstr, ".haps.allele2frequency.tsv")
  before<-read.table(prefile, sep="\t", header=F, row.names=2)
  after<-read.table(postfile, sep="\t", header=F, row.names=2)
  snp<-grep("snp:", rownames(before), fixed=T)
  beforesnp=before[snp,]
  aftersnp=after[snp,]
  #png(file=paste0("STEP11_BeforeAndAfterPhasing_", chrstr, ".png"), width=2000, height=2000, res=300)
  plot(beforesnp$V7, aftersnp$V7, xlab="Allele 2 frequency before phasing", ylab="Allele 2 frequency after phasing", main=paste0("chr ", chr))
  #dev.off()
}
dev.off()

