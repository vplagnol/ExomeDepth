
################################################################################################################################################ low level function for everted reads

countBam.everted <- function(bam.file, granges, index = bam.file, min.mapq = 1) {

  rds.counts <- numeric(length(granges))
  message('Parsing ', bam.file, ' with index ', index)
  
  rds <- Rsamtools::scanBam(file = bam.file,
                 index = index, 
                 param =Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isPaired = TRUE, isProperPair = FALSE, isNotPrimaryRead = FALSE), what = c("rname", "strand", "isize", "mapq", "pos", "isize")))[[1]]

  mapq.test <- (!is.na(rds$isize)) & (rds$mapq >= min.mapq) & !is.na(rds$pos) & (abs(rds$isize) < 100000) & ( ((rds$strand == "+") & (rds$isize < 0) ) | ((rds$strand == "-") & (rds$isize > 0) ) )
  mapq.test <- mapq.test[  !is.na(mapq.test) ]

  if (sum(mapq.test, na.rm = TRUE) > 0) {
    empty <- FALSE
    
    reads.ranges <- GenomicRanges::GRanges ( seqnames = rds$rname[ mapq.test],
                            IRanges::IRanges(start = pmin( rds$pos[ mapq.test ], rds$pos[ mapq.test ] + rds$isize [mapq.test]) , end =  pmax( rds$pos[ mapq.test ], rds$pos[ mapq.test ] + rds$isize [mapq.test])),
                            strand = rds$strand[ mapq.test ])
    
    rds.counts <- GenomicRanges::countOverlaps(granges, reads.ranges)
  }
  rds.counts      
}


################################################################################################################################################ low level function for read depth
countBamInGRanges.exomeDepth <- function (bam.file, index = bam.file, granges, min.mapq = 1, read.width = 1) {
  message("Now parsing ", bam.file)
  if (class(granges) != "GRanges") stop("Argument granges of countBamInGRanges.exomeDepth must be of the class GRanges")
  
  if (is.null(index)) index <- bam.file
  
  count.data <- rep(0, length(granges))

  
  seqs.in.target <- as.character(unique(GenomicRanges::seqnames(granges)))
  seqs.in.bam.file <- gsub(pattern = 'SN:', replacement = '',
                           as.character(grep(pattern = 'SN', unlist(Rsamtools::scanBamHeader(files = bam.file, index = index)), value = TRUE)))


  ####### first check for consistency between BAM and target regions
  if (sum(! seqs.in.target %in% seqs.in.bam.file)) {  ### if some sequences are missing
    if (sum(!paste0("chr", seqs.in.target) %in% seqs.in.bam.file) == 0) {
      warning("Apparently the BAM file uses the convention chr1 instead of 1 for chromosome names, but your target sequence does not. Therefore, adding the chr prefix to the target intervals")
      GenomicRanges::seqnames(target) <- paste0('chr', GenomicRanges::seqnames(target))
      seqs.in.target <- as.character(unique(GenomicRanges::seqnames(granges)))
    }
  }

#####  second check for consistency between BAM and target regions
  if (sum(! seqs.in.target %in% seqs.in.bam.file)) {  ### if some sequences are missing
    print("Problematic sequences:")
    print( seqs.in.target [ ! seqs.in.target %in% seqs.in.bam.file ])
    stop("Some sequences in the target data frame cannot be found in the index of the BAM file")
  }


  for (seq in seqs.in.target) { ##splits by chromosome to limit memory requirements
    #seq <- '22' ##useful for debugging
    message("Parsing chromosome ", seq)
    target.local1 <- GenomicRanges::GRanges(seqnames = seq,
                                           IRanges::IRanges(start=1, end = 5*10^8))
    
    target.local2 <-  granges[ GenomicRanges::seqnames(granges) == seq,]
    my.rows <- which (as.logical(GenomicRanges::seqnames(granges) == seq))
    
######## paired end reads
    my.param.pairs <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE, hasUnmappedMate = FALSE,
                                                isPaired = TRUE, isProperPair = TRUE, isNotPrimaryRead = FALSE),
                                              what = c("mapq", "pos", "isize"), which = target.local1)
    gal <- GenomicAlignments::readGAlignments(file = bam.file, index = index, param = my.param.pairs)
    if (length(gal) > 0) {
      gal <- as(gal, 'data.frame')
      gal <- gal[ gal$mapq > min.mapq & gal$isize > 0, ]
      
      gal <- GenomicRanges::GRanges( seqnames = gal$seqnames,
                                    IRanges::IRanges(start= gal$start, end = gal$start + gal$isize))
      count.data[ my.rows ]  <- count.data[ my.rows ] +  GenomicRanges::countOverlaps ( query = target.local2, subjec = gal)
    }

    ########### single end reads
    my.param.single <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isPaired = FALSE,  isNotPrimaryRead = FALSE),
                                               what = c("mapq", "pos"), which = target.local1)
    gal.single <- GenomicAlignments::readGAlignments(file = bam.file, index = index, param = my.param.single)
    if (length(gal.single) > 0) {
      message('Some single end reads detected in this BAM file')
      gal.single <- as(gal.single, 'data.frame')
      gal.single <- gal.single[ gal.single$mapq > min.mapq,  ]
      
      gal.single <- GenomicRanges::GRanges( seqnames = gal.single$seqnames,
                                    IRanges::IRanges(start= gal.single$start, end = gal.single$start + read.width))
      count.data[ my.rows ]  <- count.data[ my.rows ] +  GenomicRanges::countOverlaps ( query = target.local2, subjec = gal.single)
    }
  }

  if (sum(is.na(count.data) > 0)) stop('There is a bug here and count data should not contain any missing value')
  return(count.data)
}




##########################################################################################################################################################  master function for exomeDepth read counting

getBamCounts <- function(bed.frame = NULL, bed.file = NULL, bam.files, index.files = bam.files,
                         min.mapq = 20, read.width = 300, include.chr = FALSE, referenceFasta = NULL) {

  if (is.null(bed.frame)) {
    if (is.null(bed.file)) {
      stop("If no bed data frame is provided there must be a link to a bed file")
    }
    bed.frame <- read.delim(file = bed.file, header =  FALSE, stringsAsFactors = FALSE)
  }

  names(bed.frame)[1] <- 'seqnames'
  names(bed.frame)[2] <- 'start'
  names(bed.frame)[3] <- 'end'

  if (include.chr) {
    if (sum(grepl(pattern = '^chr', bed.frame$seqnames) > 0)) {
      warning('The option include.chr == TRUE adds the chr prefix to the chromosome name but it looks like the chromosome names already have a chr prefix. The argument to getBamCounts is probably an error.')
    }
    bed.frame$seqnames <- paste('chr', bed.frame$seqnames, sep = '')
  }
  
  chr.names.used <- unique(as.character(bed.frame$seqnames))
  chr.levels <- c(as.character(seq(1, 22)), subset( chr.names.used, ! chr.names.used %in% as.character(seq(1, 22))))

  bed.frame$seqnames <- factor(bed.frame$seqnames, levels = chr.levels)  ####specifying the levels is important here to not mess up the order
  bed.frame <- bed.frame[ order(bed.frame$seqnames, bed.frame$start + bed.frame$end), ]  ##order the data frame by position

  target <- GenomicRanges::GRanges(seqnames = bed.frame$seqnames,  
                    IRanges::IRanges(start=bed.frame$start+1,end=bed.frame$end))
  
  rdata <- IRanges::RangedData(space= GenomicRanges::seqnames(target),
                               ranges=GenomicRanges::ranges(target))
  
  if  ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {    
    row.names(rdata) <- make.unique(as.character(bed.frame[,4]))  ##add exon names if available
  }
  
############################################################################# add GC content
if (!is.null(referenceFasta)) {
  message('Reference fasta file provided so ExomeDepth will compute the GC content in each window')
    target.dnastringset <- Rsamtools::scanFa(referenceFasta, target)
  
    getGCcontent <- function(x) {
      GC.count <- Biostrings::letterFrequency(x,"GC")
      all.count <- Biostrings::letterFrequency(x,"ATGC")
      as.vector(ifelse(all.count==0,NA,GC.count/all.count))
    }
    rdata[["GC"]] <- getGCcontent(target.dnastringset)
  }

############################################################################# Parse BAM files
  nfiles <- length(bam.files)
  message('Parse ', nfiles, ' BAM files')
  print(bam.files)


  my.param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isPaired = TRUE, isProperPair = TRUE, isNotPrimaryRead = FALSE),
                                      what = c("mapq", "pos", "isize"), )


  
  for (i in 1:nfiles) {
    bam <- bam.files[ i ]
    index <- index.files[ i ]
    rdata[[ basename(bam) ]] <- countBamInGRanges.exomeDepth ( bam.file = bam, index = index, granges = target, min.mapq = min.mapq, read.width = read.width)
    message("Number of counted fragments : ", sum(rdata[[ basename(bam) ]]))
  }
  
  return(rdata)
}





##########################################################################################################################################################  master function for exomeDepth everted.count
count.everted.reads <- function(bed.frame = NULL,
                                bed.file = NULL,
                                bam.files,
                                index.files = bam.files,
                                min.mapq = 20,
                                include.chr = FALSE) {

  if (is.null(bed.frame)) {
    if (is.null(bed.file)) {
      stop("If no bed data frame is provided there must be a link to a bed file")
    }
    bed.frame <- read.delim(file = bed.file, header =  FALSE, stringsAsFactors = FALSE)
  }

  names(bed.frame)[1] <- 'seqnames'
  names(bed.frame)[2] <- 'start'
  names(bed.frame)[3] <- 'end'

  
  if (include.chr) bed.frame$seqnames <- paste('chr', bed.frame$seqnames, sep = '')

  target <- GenomicRanges::GRanges(seqnames = bed.frame$seqnames,  
                    IRanges::IRanges(start=bed.frame$start+1,end=bed.frame$end))
  rdata <- IRanges::RangedData(space=GenomicRanges::seqnames(target),
                               ranges=GenomicRanges::ranges(target))
  
  if  ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {    
    row.names(rdata) <- make.unique(as.character(bed.frame[,4]))  ##add exon names if available
  }

  nfiles <- length(bam.files)
  for (i in 1:nfiles) {
    bam <- bam.files[ i ]
    index <- index.files[ i ]

    rdata[[ basename(bam) ]] <- countBam.everted (bam.file = bam,
                                                  index = index,
                                                  granges = target,
                                                  min.mapq = min.mapq)
  }

  rdata <- as.data.frame(rdata)
  names(rdata)[[1]] <- 'chromosome'
  return (rdata)
}
