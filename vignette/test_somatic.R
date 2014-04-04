library(ExomeDepth)

data(exons.hg19)
data(ExomeCount)
#source('R/class_definition.R')

ExomeCount.dafr <- as(ExomeCount[,  colnames(ExomeCount)], 'data.frame')
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$space),
                                   pattern = 'chr',
                                   replacement = '')  ##remove the annoying chr letters

myTest <- somatic.CNV.call (normal  = ExomeCount.dafr$Exome1,
                            tumor = ExomeCount.dafr$Exome2,
                            prop.tumor = 0.5,
                            chromosome = ExomeCount.dafr$space,
                            start = ExomeCount.dafr$start,
                            end = ExomeCount.dafr$end,
                            names = ExomeCount.dafr$names)
