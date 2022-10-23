../samtools-1.16.1/samtools view -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00138/alignment/HG00138.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 1:2612742-2629815 1:25630000-25650000  > minimum_1_25630000_25650000.bam

../samtools-1.16.1/samtools index minimum_1_25630000_25650000.bam

mv minimum_1_25630000_25650000.bam minimum_1_25630000_25650000.bam.bai data/

#library(ExomeDepth); data(exons.hg19); my.counts <- getBamCounts(bed.frame = exons.hg19, bam.files = 'minimum_1_25630000_25650000_HG00138.bam')
#subset(my.counts, grepl(pattern = '^RHD', my.counts[['exon']]))
