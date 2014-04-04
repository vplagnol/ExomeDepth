
library(ExomeDepth)
#source('R/countBamInGranges.R')


data(exons.hg19)


my.bam <- list.files("~/Projects/exome_sequencing/IoN/Una/batch1_Una_FPD/aligned/", pattern = "*bam$", recursive = TRUE)
fasta <- "~/vyp/vincent/data/reference_genomes/fasta/human_g1k_v37.fasta"

my.counts <- getBamCounts(bed.frame = exons.hg19,
                          bam.files = my.bam,
                          referenceFasta = fasta)

