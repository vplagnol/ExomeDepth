
################################################################################################################################################ low level function for everted reads

countBam.everted <- function(bam.file, granges, index = bam.file, min.mapq = 1) {

  rds.counts <- numeric(length(granges))
  seq.names <- seqlevels(granges)
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)

  message('Parsing ', bam.file, ' with index ', index)
  
  rds <- scanBam(file = bam.file,
                 index = index, 
                 param = ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isPaired = TRUE, isProperPair = FALSE), what = c("rname", "strand", "isize", "mapq", "pos", "isize")))[[1]]

  mapq.test <- (!is.na(rds$isize)) & (rds$mapq >= min.mapq) & !is.na(rds$pos) & (abs(rds$isize) < 100000) & ( ((rds$strand == "+") & (rds$isize < 0) ) | ((rds$strand == "-") & (rds$isize > 0) ) )
  mapq.test <- mapq.test[  !is.na(mapq.test) ]

  if (sum(mapq.test) > 0) {
    empty <- FALSE
    
    reads.ranges <- GRanges ( seqnames = rds$rname[ mapq.test],
                            IRanges(start = pmin( rds$pos[ mapq.test ], rds$pos[ mapq.test ] + rds$isize [mapq.test]) , end =  pmax( rds$pos[ mapq.test ], rds$pos[ mapq.test ] + rds$isize [mapq.test])),
                            strand = rds$strand[ mapq.test ])
    
    rds.counts <- countOverlaps(granges, reads.ranges)
  }
  rds.counts      
}


################################################################################################################################################ low level function for read depth
countBamInGRanges.exomeDepth <- function (bam.file, index = bam.file, granges, min.mapq = 1, read.width = 1, force.single.end = FALSE)  {

  rds.counts <- numeric(length(granges))
  seq.names <- seqlevels(granges)
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)

  message('Parsing ', bam.file, ' with index ', index)

  for (seq.name in seq.names) {  ##for each chromosome
    if (seq.name %in% seq.names.in.bam) {
      granges.subset <- granges[seqnames(granges) == seq.name]
      strand(granges.subset) <- "*"

      empty <- TRUE
      if (!force.single.end) {
      
############################################################################# read paired end
        rds <- scanBam(file = bam.file,
                       index = index, 
                       param = ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isPaired = TRUE, isProperPair = TRUE), what = c("mapq", "pos", "isize"), which = range(granges.subset)))
        mapq.test <- (rds[[1]]$mapq >= min.mapq) & !is.na(rds[[1]]$pos) & (abs(rds[[1]]$isize) < 1000) & (rds[[1]]$isize > 0)
                                        #message('----------------- ', seq.name, ' ', length(rds[[1]]$mapq))
        
        if (sum(mapq.test) > 0 && !is.na(sum(mapq.test) )) {
          empty <- FALSE
          rds.ranges <- GRanges(seq.name, IRanges(start = rds[[1]]$pos[mapq.test], width  = rds[[1]]$isize[mapq.test]))
          rds.counts.seq.name <- countOverlaps(granges.subset, rds.ranges)
          rds.counts[as.logical(seqnames(granges) == seq.name)] <- rds.counts.seq.name
        }

      ############################################################################# read single end 
        rds <- scanBam(bam.file,
                       index = index, 
                       param = ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isPaired = FALSE), what = c("pos", "mapq", "qwidth"), which = range(granges.subset)))
        mapq.test <- (rds[[1]]$mapq >= min.mapq) & !is.na(rds[[1]]$pos) 
        
        if (sum(mapq.test) > 0 && !is.na(sum(mapq.test))) {
          empty <- FALSE
          rds.ranges <- GRanges(seq.name, IRanges(start = rds[[1]]$pos[mapq.test] - 0.5*read.width + 0.5*rds[[1]]$qwidth[ mapq.test ], width = read.width))
          rds.counts.seq.name <- countOverlaps(granges.subset, rds.ranges)
          rds.counts[as.logical(seqnames(granges) == seq.name)] <- rds.counts.seq.name
        }
      }

      if (force.single.end) {  ##request to deal with reads in a single end manner
        rds <- scanBam(bam.file,
                       index = index,
                       param = ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE), what = c("pos", "mapq", "qwidth"), which = range(granges.subset)))

        mapq.test <- (rds[[1]]$mapq >= min.mapq) & !is.na(rds[[1]]$pos)
        if (sum(mapq.test) > 0 && !is.na(sum(mapq.test))) {
          empty <- FALSE
          rds.ranges <- GRanges(seq.name, IRanges(start = rds[[1]]$pos[mapq.test] - 0.5*read.width + 0.5*rds[[1]]$qwidth[ mapq.test ], width = read.width))
          rds.counts.seq.name <- countOverlaps(granges.subset, rds.ranges)
          rds.counts[as.logical(seqnames(granges) == seq.name)] <- rds.counts.seq.name
        }
      }
      
######################
      if (empty) rds.counts[as.logical(seqnames(granges) == seq.name)] <- 0  ## do I need that?
                                        #message('Sequence ', seq.name, ' ',  rds.counts[as.logical(seqnames(granges) == seq.name)], '\n')
    }
  }
  rds.counts
}





##########################################################################################################################################################  master function for exomeDepth read counting

getBamCounts <- function(bed.frame = NULL, bed.file = NULL, bam.files, index.files = bam.files,
                         min.mapq = 20, read.width = 300, include.chr = FALSE, referenceFasta = NULL, force.single.end = FALSE) {
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
  chr.names.used <- unique(as.character(bed.frame$seqnames))
  chr.levels <- c(as.character(seq(1, 22)), subset( chr.names.used, ! chr.names.used %in% as.character(seq(1, 22))))

  bed.frame$seqnames <- factor(bed.frame$seqnames, levels = chr.levels)  ####specifying the levels is important here to not mess up the order
  bed.frame <- bed.frame[ order(bed.frame$seqnames, bed.frame$start + bed.frame$end), ]  ##order the data frame by position

  target <- GRanges(seqnames = bed.frame$seqnames,  
                    IRanges(start=bed.frame$start+1,end=bed.frame$end))
  
  rdata <- RangedData(space=seqnames(target),
                      ranges=ranges(target))
  
  if  ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {    
    row.names(rdata) <- make.unique(as.character(bed.frame[,4]))  ##add exon names if available
  }
  
############################################################################# add GC content
if (!is.null(referenceFasta)) {
  message('Reference fasta file provided so exomeDepth will compute the GC content in each window')
    target.dnastringset <- scanFa(referenceFasta, target)
  
    getGCcontent <- function(x) {
      GC.count <- letterFrequency(x,"GC")
      all.count <- letterFrequency(x,"ATGC")
      as.vector(ifelse(all.count==0,NA,GC.count/all.count))
    }
    rdata[["GC"]] <- getGCcontent(target.dnastringset)
  }

############################################################################# Parse BAM files
  nfiles <- length(bam.files)
  message('Parse ', nfiles, ' BAM files')
  print(bam.files)
  
  for (i in 1:nfiles) {
    bam <- bam.files[ i ]
    index <- index.files[ i ]
    
    rdata[[ basename(bam) ]] <- countBamInGRanges.exomeDepth(bam.file = bam,
                                                             index = index,
                                                             granges = target,
                                                             min.mapq = min.mapq,
                                                             read.width = read.width,
                                                             force.single.end = force.single.end)
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

  target <- GRanges(seqnames = bed.frame$seqnames,  
                    IRanges(start=bed.frame$start+1,end=bed.frame$end))
  rdata <- RangedData(space=seqnames(target),
                      ranges=ranges(target))
  
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
