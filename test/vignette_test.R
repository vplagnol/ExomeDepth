library(GenomicRanges)
library(VGAM)
library(aod)

library(ExomeDepth)
source("R/class_definition.R")
source("R/tools.R")
source("R/optimize_reference_set.R")
dyn.load('src/ExomeDepth.so')




load('data/ExomeCount.RData')
data <- as(ExomeCount, 'data.frame')
data$space[1:10000] <- 2

test <- new('ExomeDepth',
            test = data$Exome2,
            reference = data$Exome4,
            formula = 'cbind(test, reference) ~ 1',
            subset.for.speed = seq(1, nrow(data), 100))

show(test)



my.test <- data$Exome3
my.reference.set <- as.matrix(data[, c('Exome1', 'Exome2', 'Exome4')])
my.choice <- select.reference.set (test.counts = my.test, reference.counts = my.reference.set, bin.length = (data$end - data$start)/1000, n.bins.reduced = 10000)
#print(my.choice)

my.reference.selected <- apply(as.matrix( data[, my.choice$reference.choice] ), MAR = 1, FUN = sum)


test <- new('ExomeDepth',
            test = my.test,
            reference = my.reference.selected,
            formula = 'cbind(test, reference) ~ 1')

#source('R/call_CNVs_method.R')
my.calls <- CallCNVs(test, transition.probability = 10^-4, chromosome = data$space, start = data$start, end = data$end, name = data$names)
print(my.calls@CNV.calls)
