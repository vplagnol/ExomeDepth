\name{countBamInGRanges.exomeDepth}
\alias{countBamInGRanges.exomeDepth}
\title{
  Compute read count data from BAM files.
}
\description{
  Parses a BAM file and count reads that are located within a target
  region defined by a GenomicRanges object.
}
\usage{
countBamInGRanges.exomeDepth(bam.file, index = bam.file, granges,
min.mapq = 1, read.width = 1, force.single.end = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam.file}{BAM file to be parsed}
  \item{index}{Index of the BAM file, without the '.bai'
    suffix.}
  \item{granges}{Genomic ranges object defining the bins}
  \item{min.mapq}{Minimum read mapping quality (Phred scaled).}
  \item{read.width}{For single end reads, an estimate of the frament
    size. For paired reads, the fragment size can be directly computed
    from the paired alignment and this value is ignored.}
  \item{force.single.end}{Treat all paired reads as single end. It is
    generally a bad idea but may be useful in some rare situations. Only
    use this option if you really know what you are doing.}
}
\details{
  Largely derived from its equivalent function in the exomeCopy package.
}
\value{
  A GRanges object with count data.
}
