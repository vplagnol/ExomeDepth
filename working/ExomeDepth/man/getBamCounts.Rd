\name{getBamCounts}
\alias{getBamCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get count data for multiple exomes
}
\description{
Essentially a wrapper for the accessory function countBamInGRanges which only considers a single BAM file at a time. 
}
\usage{
getBamCounts(bed.frame = NULL,
bed.file = NULL,
bam.files,
index.files = bam.files,
min.mapq = 20,
read.width = 300,
include.chr = FALSE,
referenceFasta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed.frame}{\code{data.frame} containing the definition of the
    regions.
    The first three columns must be chromosome, start, end.
}
\item{bed.file}{\code{character} file name. Target BED file with the
  definition of the regions. This file will only be used if no bed.frame
  argument is provided. No headers are assumed so remove them if they
  exist. Either a bed.file or a bed.frame must be provided for this
  function to run.}
  \item{bam.files}{
    \code{character}, list of BAM files to extract read count data from.
  }
  \item{index.files}{
    Optional \code{character} argument with the list of indexes for the BAM files, without the
    '.bai' suffix. If the indexes are simply obtained by adding .bai to
    the BAM files, this argument does not need to be specified.
  }
\item{min.mapq}{
    \code{numeric}, minimum mapping quality to include a read.
}
  \item{read.width}{
    \code{numeric}, maximum distance between the side of the target region and the middle of the paired read
    to include the paired read into that region.
}
  \item{include.chr}{
    \code{logical}, if set to TRUE, this function will remove the string 'chr' from the
    chromosome names of the target BED file.
  }
  \item{referenceFasta}{
    \code{character}, file name for the reference genome in fasta
    format. If available, GC content will be computed and added to the
    output.}
}
\details{
  This function is largely a copy of a similar one available in the
  exomeCopy package.
}
\value{
  A GenomicRanges object that stores the read count data for the BAM
  files listed as argument.
}
\references{
  exomeCopy R package.
}
\author{
  Vincent Plagnol
}
\examples{
\dontrun{
load(exons.hg19)

my.counts <- getBamCounts(bed.frame = exonpos,
                          bam.files = my.bam,
                          referenceFasta = 'human_g1k_v37.fasta')  
} 
}