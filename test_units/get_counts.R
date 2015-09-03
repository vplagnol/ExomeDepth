library(ExomeDepth)

my.bam <- "/scratch2/vyp-scratch2/exomes_temp/Shamima/Shamima_April2014/aligned/SRBEI_16/SRBEI_16_sorted_unique.bam"
fasta <- "/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta"

data(exons.hg19)

bed.frame <- exons.hg19
names(bed.frame)[1] <- 'seqnames'
names(bed.frame)[2] <- 'start'
names(bed.frame)[3] <- 'end'

source('R/countBamInGranges.R')
granges <- GenomicRanges::GRanges(seqnames = bed.frame$seqnames,
                                  IRanges::IRanges(start=bed.frame$start+1,end=bed.frame$end))

test <- countBamInGRanges.exomeDepth (bam.file = my.bam, granges = granges, min.mapq = 10)
save(list = 'test', file = 'overlap2.RData')

