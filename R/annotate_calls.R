

overlap.feature <- function(chromosome, start, end, index.chromosome, index.start, index.end, index.name, fraction.overlap.required) {

  
  return( .Call("C_overlap", chromosome, start, end, index.chromosome, index.start, index.end, index.name, fraction.overlap.required))

}


test <- function() {

  dyn.load('src/overlap.so')
  Conrad <- read.table('/ugi/home/shared/vincent/annotations/Conrad_CNVs_hg19.tab', header = TRUE, stringsAsFactors = FALSE)
  my.res <- overlap.feature(chromosome = '1', start = 10^5, end = 3*10^5, index.chromosome = Conrad$chr, index.start = Conrad$start, index.end = Conrad$end, index.name = Conrad$CNVR, fraction.overlap.required = 0.5)

}
test()
