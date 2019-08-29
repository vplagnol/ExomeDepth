library(GenomicRanges)


#############################################################################
#load('../../CNVcalls/depth_Sergey_AROS.RData')
#ExomeCount <- my.counts[space(my.counts) == '1',]
#colnames(ExomeCount) <- c('GC', 'Exome1', 'Exome2', 'Exome3', 'Exome4')
#save(list = 'ExomeCount', file = 'data/ExomeCount.RData')


system("sort -k1,1 -k2,2n /cluster/project8/vyp/vincent/toolsVarious/ensemblAPI/data/canonical_all_genes_Human_hg19.bed | ~/Software/bedtools-2.17.0/bin/bedtools merge  -d 100 -i - -nms   | sed -e 's/;/,/g'  > data/bedFiles/genes_hg19.bed")
system("sort -k1,1 -k2,2n /cluster/project8/vyp/vincent/toolsVarious/ensemblAPI/data/canonical_all_exons_Human_hg19.bed | ~/Software/bedtools-2.17.0/bin/bedtools merge  -d 50 -i - -nms   | sed -e 's/;/,/g'  > data/bedFiles/exons_hg19.bed" )



genes.hg19 <- read.table('data/bedFiles/genes_hg19.bed', col.names = c('chromosome', 'start', 'end', 'name'), stringsAsFactors = FALSE)
save(list = c('genes.hg19') , file = 'data/genes.hg19.RData')




#exons.hg19 <- read.table('~/vyp/vincent/data/reference_genomes/annotations/exonPos_hg19.tab', col.names = c('chromosome', 'start', 'end', 'name'))


exons.hg19 <- read.table('data/bedFiles/exons_hg19.bed', col.names = c('chromosome', 'start', 'end', 'name'), stringsAsFactors = FALSE)    
exons.hg19.X <- subset(exons.hg19, chromosome %in% c('X'))
exons.hg19.Y <- subset(exons.hg19, chromosome %in% c('Y'))

exons.hg19 <- subset(exons.hg19, !(chromosome %in% c('X', 'Y')))
exons.hg19$name <- as.character(exons.hg19$name)
save(list = c('exons.hg19') , file = 'data/exons.hg19.RData')
save(list = c('exons.hg19.X') , file = 'data/exons.hg19.X.RData')


data <- read.table('/cluster/project8/vyp/vincent/data/reference_genomes/annotations/Conrad_CNVs_hg19.tab', header = TRUE)
Conrad.hg19.common.CNVs <- GRanges(seqnames = data$chr,
                                   IRanges(start=data$start,end=data$end),
                                   names = data$CNVR)
save(list = 'Conrad.hg19.common.CNVs', file = 'data/Conrad.hg19.RData')






