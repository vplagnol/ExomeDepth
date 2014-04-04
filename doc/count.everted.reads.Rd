\name{count.everted.reads}
\alias{count.everted.reads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Count the number of everted reads for a set of BAM files.
}
\description{
  This is the ExomeDepth high level function that takes a GenomicRanges
  object, a list of indexed/sorted BAM files, and compute the number of
  everted reads in each of the defined bins.
}
\usage{
count.everted.reads(bed.frame = NULL, bed.file = NULL,
                    bam.files, index.files = bam.files,
                    min.mapq = 20, include.chr = FALSE)
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
  \item{include.chr}{
    \code{logical}, if set to TRUE, this function will remove the string 'chr' from the
    chromosome names of the target BED file.
  }
}
\value{
  A data frame that contains the region and the number of identified
  reads in each bin.
}
\details{
  Everted reads are characteristic of the presence of duplications in a
  BAM files. This routine will parse a BAM files and the suggested use
  is to provide relatively large bins (for example gene based, and
  ExomeDepth has a genes.hg19 object that is appropriate for this) to
  flag the genes that contain such reads suggestive of a duplication. A
  manual check of the data using IGV is recommended to confirm that
  these reads are all located in the same DNA region, which would
  confirm the presence of a copy number variant.
}
\references{
  Computational methods for discovering structural variation with
  next-generation sequencing, Medvedev P, Stanciu M, Brudno M., Nature
  Methods 2009
}
\author{
  Vincent Plagnol
}
\note{
  This function calls a lower level function called XXX that works on each
single BAM file.
}
\seealso{
  getBAMCounts
}
\examples{
\dontrun{  test <- count.everted.reads (bed.frame = genes.hg19,
  bed.file = NULL,
  bam.files = bam.files,
  min.mapq = 20,
  include.chr = FALSE)
}
}

