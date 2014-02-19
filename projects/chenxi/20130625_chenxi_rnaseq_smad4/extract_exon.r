library(GenomicFeatures)
library(DEXSeq)
hse <- makeTranscriptDbFromBiomart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
exonicParts <- prepareAnnotationForDEXSeq( hse )
bamDir <- system.file("extdata",package="parathyroidSE",mustWork=TRUE)
fls <- list.files(bamDir, pattern="bam$",full=TRUE)
bamlst <- BamFileList(fls)
SE <- countReadsForDEXSeq( exonicParts, bamlst )
ecs <- buildExonCountSet( SE, c("A", "A", "B"), exonicParts )