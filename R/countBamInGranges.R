#' Positions of exons on build hg19 of the human genome
#'
#' Exon position extracted from the ensembl database version 71.
#'
#' @name exons.hg19
#' @md
#' @docType data
#' @format A data frame with 192,379 observations on the following 4 variables:
#' * chromosome, a factor with levels `1`, `2` `3` `4`, `5` `6` `7` `8` `9`, `10` `11` `12` `13` `14` `15` `16` `17` `18` `19` `2` `20` `21` `22`
#' * start a numeric vector
#' * end a numeric vector
#' * name A character vector of names for the exon(s)
#' @source Ensemble database version 71.
#' @keywords datasets
NULL

#' Positions of exons on build hg19 of the human genome and on chromosome X
#'
#' Exon position extracted from the ensembl database version 61 and on
#' chromosome X only.
#'
#' @name exons.hg19.X
#' @md
#' @docType data
#' @format A data frame of exons with the following 4 variables:
#' * chromosome, a factor with levels `X`, `Y`.
#' * start Numeric.
#' * end Numeric.
#' * name Character names for the exons.
#' @source Ensemble database version 71.
#' @keywords datasets
NULL


#' Positions of genes on build hg19 of the human genome
#'
#' Exon position extracted from the ensembl database version 71.
#'
#'
#' @name genes.hg19
#' @docType data
#' @md
#' @format A data frame with 18,033 observations on the following 4 variables:
#' * chromosome, a factor with levels `1`, `2` `3` `4`, `5` `6` `7` `8` `9`, `10` `11` `12` `13` `14` `15` `16` `17` `18` `19` `2` `20` `21` `22`
#' * start a numeric vector
#' * end a numeric vector
#' * name A character vector of names for the exon(s)
#' @source Ensemble database version 71.
#' @keywords datasets
NULL


#' Conrad et al common CNVs
#'
#' Positions of common CNV calls (detected in a panel of 42 sample) from the
#' Conrad et al paper (Nature 2010). This is build hg19 of the human genome.
#'
#'
#' @name Conrad.hg19.common.CNVs
#' @docType data
#' @format A data frame with common CNV calls.
#' @source Conrad et al, Origins and functional impact of copy number variation in the human genome, Nature 2010
#' @keywords datasets
NULL


#' Example dataset for ExomeDepth
#'
#' An example dataset of 4 exome samples, chromosome 1 only.
#'
#'
#' @name ExomeCount
#' @docType data
#' @md
#' @format A data frame with 25592 observations on the following 9 variables:
#' * chromosome, Character vector with chromosome names (only chromosome 1 in that case)
#' * start, start of exons
#' * end, end of exons
#' * exons, character name of exons
#' * camfid.032KA_sorted_unique.bam
#' * camfid.033ahw_sorted_unique.bam
#' * camfid.034pc_sorted_unique.bam
#' * camfid.035if_sorted_unique.bam
#' * GC, a numeric vector with the GC content
#' @source Dataset generated in collaboration with Sergey Nejentsev, University of Cambridge.
#' @keywords datasets
NULL


#' Counts everted reads from a single BAM file
#'
#' This is a utility function that is called by the higher level
#' count.everted.reads. It processes each BAM file individually to generate the
#' count data.
#'
#' Most users will not use this function, and it will only be called by the
#' higher level count.everted.reads. Nevertheless it may be useful on its own
#' in some cases.
#'
#' @param bam.file BAM file that needs to be parsed
#' @param granges Genomic Ranges object with the location of the bins for which
#' we want to count the everted reads.
#' @param index Index for the BAM files.
#' @param min.mapq Minimum mapping quality to include reads.
#' @return A list with the number of reads in each bin.
#' @seealso count.everted.reads


countBam.everted <- function(bam.file, granges, index = bam.file, min.mapq = 1) {

  rds.counts <- numeric(length(granges))
  message('Parsing ', bam.file, ' with index ', index)

  rds <- Rsamtools::scanBam(file = bam.file,
                 index = index,
                 param =Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isPaired = TRUE, isProperPair = FALSE, isSecondaryAlignment = FALSE), what = c("rname", "strand", "isize", "mapq", "pos", "isize")))[[1]]

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



#' Compute read count data from BAM files.
#'
#' Parses a BAM file and count reads that are located within a target region
#' defined by a GenomicRanges object.
#'
#' Largely derived from its equivalent function in the exomeCopy package.
#'
#' @param bam.file BAM file to be parsed
#' @param index Index of the BAM file, without the '.bai' suffix.
#' @param granges Genomic ranges object defining the bins
#' @param min.mapq Minimum read mapping quality (Phred scaled).
#' @param read.width For single end reads, an estimate of the frament size. For
#' paired reads, the fragment size can be directly computed from the paired
#' alignment and this value is ignored.
#' @return A GRanges object with count data.

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
    print(seqs.in.bam.file)
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
                                                isPaired = TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE),
                                              what = c("mapq", "pos", "isize"), which = target.local1)
    gal <- GenomicAlignments::readGAlignments(file = bam.file, index = index, param = my.param.pairs)
    if (length(gal) > 0) {
      gal <- methods::as(gal, 'data.frame')
      gal <- gal[ gal$mapq > min.mapq & gal$isize > 0, ]

      gal <- GenomicRanges::GRanges( seqnames = gal$seqnames,
                                    IRanges::IRanges(start= gal$start, end = gal$start + gal$isize))
      count.data[ my.rows ]  <- count.data[ my.rows ] +  GenomicRanges::countOverlaps ( query = target.local2, subjec = gal)
    }

    ########### single end reads
    my.param.single <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isPaired = FALSE,  isSecondaryAlignment = FALSE),
                                               what = c("mapq", "pos"), which = target.local1)
    gal.single <- GenomicAlignments::readGAlignments(file = bam.file, index = index, param = my.param.single)
    if (length(gal.single) > 0) {
      message('Some single end reads detected in this BAM file')
      gal.single <- methods::as(gal.single, 'data.frame')
      gal.single <- gal.single[ gal.single$mapq > min.mapq,  ]

      gal.single <- GenomicRanges::GRanges( seqnames = gal.single$seqnames,
                                    IRanges::IRanges(start= gal.single$start, end = gal.single$start + read.width))
      count.data[ my.rows ]  <- count.data[ my.rows ] +  GenomicRanges::countOverlaps ( query = target.local2, subjec = gal.single)
    }
  }

  if (sum(is.na(count.data) > 0)) stop('There is a bug here and count data should not contain any missing value')
  return(count.data)
}





#' Get count data for multiple exomes
#'
#' Essentially a wrapper for the accessory function countBamInGRanges which
#' only considers a single BAM file at a time.
#'
#' This function is largely a copy of a similar one available in the exomeCopy
#' package.
#'
#' @param bed.frame \code{data.frame} containing the definition of the regions.
#' The first three columns must be chromosome, start, end.
#' @param bed.file \code{character} file name. Target BED file with the
#' definition of the regions. This file will only be used if no bed.frame
#' argument is provided. No headers are assumed so remove them if they exist.
#' Either a bed.file or a bed.frame must be provided for this function to run.
#' @param bam.files \code{character}, list of BAM files to extract read count
#' data from.
#' @param index.files Optional \code{character} argument with the list of
#' indexes for the BAM files, without the '.bai' suffix. If the indexes are
#' simply obtained by adding .bai to the BAM files, this argument does not need
#' to be specified.
#' @param min.mapq \code{numeric}, minimum mapping quality to include a read.
#' @param read.width \code{numeric}, maximum distance between the side of the
#' target region and the middle of the paired read to include the paired read
#' into that region.
#' @param include.chr \code{logical}, if set to TRUE, this function will add
#' the string 'chr' to the chromosome names of the target BED file.
#' @param referenceFasta \code{character}, file name for the reference genome
#' in fasta format. If available, GC content will be computed and added to the
#' output.
#' @return A GenomicRanges object that stores the read count data for the BAM
#' files listed as argument.
#' @author Vincent Plagnol
#' @references exomeCopy R package.
#' @examples
#'
#' \dontrun{
#' load(exons.hg19)
#'
#' my.counts <- getBamCounts(bed.frame = exonpos,
#'                           bam.files = my.bam,
#'                           referenceFasta = 'human_g1k_v37.fasta')
#' }
#'

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

  if  ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {
      GenomicRanges::values(target) <- cbind(GenomicRanges::values(target), data.frame(exon = as.character(bed.frame[,4]),stringsAsFactors = FALSE))
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
  GenomicRanges::values(target) <- cbind(GenomicRanges::values(target), data.frame(GC = getGCcontent(target.dnastringset)))
}

############################################################################# Parse BAM files
  nfiles <- length(bam.files)
  message('Parse ', nfiles, ' BAM files')
  print(bam.files)


  my.param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isPaired = TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE),
                                      what = c("mapq", "pos", "isize"), )


  exon_count_frame <- dplyr::tibble(chromosome = as(GenomicRanges::seqnames(target), 'character'),
                                    start = as(GenomicRanges::start(target), 'numeric'),
                                    end = as(GenomicRanges::end(target), 'numeric'))
  exon_count_frame <- dplyr::bind_cols(exon_count_frame, as(GenomicRanges::values(target), 'data.frame'))

  for (i in 1:nfiles) {
    bam <- bam.files[ i ]
    index <- index.files[ i ]
    exon_count_frame[[ basename(bam) ]] <- countBamInGRanges.exomeDepth ( bam.file = bam, index = index, granges = target, min.mapq = min.mapq, read.width = read.width)
    message("Number of counted fragments : ", sum(exon_count_frame[[ basename(bam) ]]))
  }


  return(data.frame(exon_count_frame))
}






#' Count the number of everted reads for a set of BAM files.
#'
#' This is the ExomeDepth high level function that takes a GenomicRanges
#' object, a list of indexed/sorted BAM files, and compute the number of
#' everted reads in each of the defined bins.
#'
#' Everted reads are characteristic of the presence of duplications in a BAM
#' files. This routine will parse a BAM files and the suggested use is to
#' provide relatively large bins (for example gene based, and ExomeDepth has a
#' genes.hg19 object that is appropriate for this) to flag the genes that
#' contain such reads suggestive of a duplication. A manual check of the data
#' using IGV is recommended to confirm that these reads are all located in the
#' same DNA region, which would confirm the presence of a copy number variant.
#'
#' @param bed.frame \code{data.frame} containing the definition of the regions.
#' The first three columns must be chromosome, start, end.
#' @param bed.file \code{character} file name. Target BED file with the
#' definition of the regions. This file will only be used if no bed.frame
#' argument is provided. No headers are assumed so remove them if they exist.
#' Either a bed.file or a bed.frame must be provided for this function to run.
#' @param bam.files \code{character}, list of BAM files to extract read count
#' data from.
#' @param index.files Optional \code{character} argument with the list of
#' indexes for the BAM files, without the '.bai' suffix. If the indexes are
#' simply obtained by adding .bai to the BAM files, this argument does not need
#' to be specified.
#' @param min.mapq \code{numeric}, minimum mapping quality to include a read.
#' @param include.chr \code{logical}, if set to TRUE, this function will add
#' the string 'chr' to the chromosome names of the target BED file.
#' @return A data frame that contains the region and the number of identified
#' reads in each bin.
#' @note This function calls a lower level function called XXX that works on
#' each single BAM file.
#' @seealso getBAMCounts
#' @references Computational methods for discovering structural variation with
#' next-generation sequencing, Medvedev P, Stanciu M, Brudno M., Nature Methods
#' 2009
#' @examples
#'
#' \dontrun{  test <- count.everted.reads (bed.frame = genes.hg19,
#'   bed.file = NULL,
#'   bam.files = bam.files,
#'   min.mapq = 20,
#'   include.chr = FALSE)
#' }
#'

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

  if  ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {
    GenomicRanges::values(target) <- cbind(GenomicRanges::values(target), data.frame(exon = as.character(bed.frame[,4]),stringsAsFactors = FALSE))
  }

  exon_count_frame <- dplyr::tibble(chromosome = as(GenomicRanges::seqnames(target), 'character'),
                                    start = as(GenomicRanges::start(target), 'numeric'),
                                    end = as(GenomicRanges::end(target), 'numeric'))
  exon_count_frame <- dplyr::bind_cols(exon_count_frame, as(GenomicRanges::values(target), 'data.frame'))

  nfiles <- length(bam.files)
  for (i in 1:nfiles) {
    bam <- bam.files[ i ]
    index <- index.files[ i ]
    exon_count_frame[[ basename(bam) ]] <- countBam.everted (bam.file = bam,  ## replace old RangedData with GRanges
                                                  index = index,
                                                  granges = target,
                                                  min.mapq = min.mapq)
   }

  return(data.frame(exon_count_frame))
}
