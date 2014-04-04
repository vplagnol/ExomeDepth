\name{countBam.everted}
\alias{countBam.everted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Counts everted reads from a single BAM file
}
\description{
  This is a utility function that is called by the higher level
  count.everted.reads. It processes each BAM file individually to
  generate the count data.
}
\usage{
countBam.everted(bam.file, granges, index = bam.file, min.mapq = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam.file}{
    BAM file that needs to be parsed
}
\item{granges}{
  Genomic Ranges object with the location of the bins for which we want
  to count the everted reads.
}
  \item{index}{
    Index for the BAM files.
}
  \item{min.mapq}{
    Minimum mapping quality to include reads.
}
}
\details{
  Most users will not use this function, and it will only be called by
  the higher level count.everted.reads. Nevertheless it may be useful on
  its own in some cases.
}
\value{
  A list with the number of reads in each bin.
}
\author{
Vincent Plagnol
}
\seealso{
  count.everted.reads
}
