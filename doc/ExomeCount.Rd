\name{ExomeCount}
\alias{ExomeCount}
\docType{data}
\title{
Example dataset for ExomeDepth
}
\description{
An example dataset of 4 exome samples, chromosome 1 only.
}
\usage{data(ExomeCount)}
\format{
  A data frame with 25592 observations on the following 9 variables.
  \describe{
    \item{\code{chromosome}}{a character vector with chromosome names (only chromosome 1 in that case)}
    \item{\code{start}}{Start of the exons}
    \item{\code{end}}{End of the exons}
    \item{\code{exons}}{a character vector with the exon names.}
    \item{\code{camfid.032KA_sorted_unique.bam}}{a numeric vector of read count data.}
    \item{\code{camfid.033ahw_sorted_unique.bam}}{a numeric vector of read count data.}
    \item{\code{camfid.035if_sorted_unique.bam}}{a numeric vector of read count data.}
    \item{\code{camfid.034pc_sorted_unique.bam}}{a numeric vector of read count data.}
    \item{\code{GC}}{a numeric vector with the GC content.}
  }
}
\source{
  Dataset generated in collaboration with Sergey Nejentsev, University
  of Cambridge.
}
\references{
  Paper currently under review.
}
\keyword{datasets}
