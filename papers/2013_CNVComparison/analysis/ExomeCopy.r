library(exomeCopy)

setwd("/scratch/cqs/shengq1/dnaseq/2110/ExomeCopy/result")

sample.names <- c(
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

my.bed<-"/scratch/cqs/lij17/cnv/SureSelect_XT_Human_All_Exon_V4_withoutchr_withoutY_lite.bed"

target.df <- read.delim(my.bed, header = TRUE)

target <- GRanges(seqname = target.df[, 1], IRanges(start = target.df[,2] , end = target.df[, 3]))

counts <- RangedData(space = seqnames(target), ranges = ranges(target))

for (i in 1:length(my.bams)) {
  counts[[sample.names[i]]] <- countBamInGRanges(my.bams[i],
                                                 target)
}

counts[["GC"]] <- getGCcontent(target, my.fasta)

counts[["GC.sq"]] <- example.counts$GC^2

counts[["bg"]] <- generateBackground(sample.names, counts, median)

counts[["width"]] <- width(example.counts)

fit.list <- lapply(sample.names, function(sample.name) {
  lapply(seqlevels(target), function(seq.name) {
    exomeCopy(counts[seq.name], sample.name, X.names = c("bg", "GC", "GC.sq", "width"), S = 0:4, d = 2)
  })
})

compiled.segments <- compileCopyCountSegments(fit.list)

save(compiled.segments, file="compiled.segments.RData")

CNV.segments <- compiled.segments[compiled.segments$copy.count != 2, ]

CNV.segments <- CNV.segments[CNV.segments$nranges > 5, ]

CNV.df<-as.data.frame(CNV.segments)

write.table(x=CNV.df, file = "ExomeCopy.tsv", row.names = FALSE, sep="\t")


