## migration code

load("bak/ExomeCount.RData")

count_values <- ExomeCount@values[[1]]
my_nrow <- nrow(count_values)

my_frame <- GenomicRanges::GRanges(rep("chr1", my_nrow),
	 ranges =  IRanges::IRanges(start = start(ExomeCount),
                   	end = end(ExomeCount)))


for (i in 1:ncol(count_values)) {
   elementMetadata(my_frame)[[ names(count_values)[i] ]] <- count_values[[i]]
}
elementMetadata(my_frame)[['names']] <- names(ExomeCount@ranges@unlistData)

ExomeCount <- my_frame


save(list = 'ExomeCount', file = 'data/ExomeCount.RData')

