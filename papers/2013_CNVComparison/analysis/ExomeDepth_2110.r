library(ExomeDepth)

setwd("/scratch/cqs/shengq1/dnaseq/2110/ExomeDepth/result")

data(exons.hg19)

SampleNames <- c(
  "2110_JP_01"
  ,"2110_JP_02"
  ,"2110_JP_03"
  ,"2110_JP_04"
  ,"2110_JP_05"
  ,"2110_JP_06"
  ,"2110_JP_07"
  ,"2110_JP_08"
  ,"2110_JP_09"
  ,"2110_JP_10"
  ,"2110_JP_11"
  ,"2110_JP_12"
  ,"2110_JP_13"
  ,"2110_JP_14"
  ,"2110_JP_15"
  ,"2110_JP_16"
)

my.bams <- c(
  "/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-1_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-2_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-3_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-4_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-5_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-6_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-7_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-8_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-9_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-10_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-11_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-12_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-13_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-14_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-15_realigned_recal_rmdup.sorted.bam"
  ,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-16_realigned_recal_rmdup.sorted.bam"
)

my.fasta<-"/data/cqs/guoy1/reference/hg19/hg19_chr.fa"

countfile<-"ExomeCount.RData"

ExomeCount <- getBamCounts(bed.frame = exons.hg19,
                           bam.files = my.bams,
                           include.chr = FALSE,
                           referenceFasta = my.fasta)

save(ExomeCount, file=countfile)

ExomeCount.dafr <- as(ExomeCount[, colnames(ExomeCount)], 'data.frame')

colnames(ExomeCount.dafr)[7:ncol(ExomeCount.dafr)]<-SampleNames

ExomeCount.mat <- as.matrix(ExomeCount.dafr[, c(7:ncol(ExomeCount.dafr))])

nsamples <- ncol(ExomeCount.mat)

### start looping over each sample
for (i in 1:nsamples) {
  
  #### Create the aggregate reference set for this sample
  my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
                                     reference.counts = ExomeCount.mat[,-i],
                                     bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                     n.bins.reduced = 10000)
  
  my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)
  
  message('Now creating the ExomeDepth object')
  all.exons <- new('ExomeDepth',
                   test = ExomeCount.mat[,i],
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')
  
  ################ Now call the CNVs
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = ExomeCount.dafr$space,
                        start = ExomeCount.dafr$start,
                        end = ExomeCount.dafr$end,
                        name = ExomeCount.dafr$names)
  
  
  output.file <- paste0(colnames(ExomeCount.mat)[i], '.tsv')
  write.table(file = output.file, x = all.exons@CNV.calls, row.names = FALSE, sep="\t")
}

