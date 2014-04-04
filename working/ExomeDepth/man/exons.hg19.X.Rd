\name{exons.hg19.X}
\alias{exons.hg19.X}
\docType{data}
\title{
  Positions of exons on build hg19 of the human genome and on chromosome X
}
\description{
  Exon position extracted from the ensembl database version 61 and on
  chromosome X only.
}
\usage{data(exons.hg19)}
\format{
  A data frame of exons with the following 4 variables.
  \describe{
    \item{\code{chromosome}}{a factor with levels \code{X} \code{Y}}
    \item{\code{start}}{a numeric vector}
    \item{\code{end}}{a numeric vector}
    \item{\code{name}}{A character vector of names for the exons}
  }
}
\source{
  Ensemble database version 71.
}
\keyword{datasets}
