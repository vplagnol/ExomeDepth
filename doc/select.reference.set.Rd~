\name{select.reference.set}
\alias{select.reference.set}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine multiple samples to optimize the reference set in order to maximise the power to detect CNV.
}
\description{
The power to detect copy number variant (CNVs) from targeted sequence data can be
maximised if the most appropriate set of sequences is used as
reference. This function is designed to combine multiple reference
exomes in order to build the best reference set.
}
\usage{
select.reference.set(test.counts, reference.counts, bin.length = NULL,
n.bins.reduced = 0, data = NULL, formula = 'cbind(test, reference) ~ 1',
phi.bins = 1)
}
\arguments{
  \item{test.counts}{
    Read count data for the test sample (numeric, typically a vector of
    integer values).
}
  \item{reference.counts}{
    Matrix of read count data for a set of additional samples that can
    be used as a comparison point for the test sample.
}
  \item{bin.length}{
    Length (in bp) of each of the regions (often exons, but not
    necessarily) that were used to compute the read count data
    (i.e. what is provided in the argument test.counts of this
    function). 
    If not provided all bins are assumed to have equal length.
  }
  \item{n.bins.reduced}{
    This function can be slow when applied genome-wide. For the purpose
    of building the reference sample, it is not necessary to use the
    full data. The number provided by this argument specifies the number
    of regions (typically exons) that will be sub-sampled (using a grid)
    to optimise the referenceset. We find that 10,000 is largely
    sufficient for exome data.
  }
  \item{data}{Defaults to NULL: A data frame of covariates that can be included in the
    model.}
  \item{formula}{Defaults to 'cbind(test, reference) ~ 1'. This formula
    will be used to fit the read count data. Covariates present in the
    data frame (for example GC content) can be included in the right
    hand side of the equation'. If covariates are provided they must be
    provided as arguments (in the data frame ``data'').}
  \item{phi.bins}{Numeric integer (typically 1, 2, or 3) that specifies the number of windows
    where the over-dispersion parameter phi can vary. It defaults to 1,
    i.e. a single over-dispersion parameter, independently of read depth.}
}
\value{
 \item{reference.choice }{character: list of samples selected as optimum
 reference set.}
 \item{summary.stats}{A data frame summarizing the output of this
   computation, including expected Bayes factor, Rs statistic (see
   reference for explanation) for multiple choices of reference set.}
}
\references{
  Key paper currently under review.
}
\author{
  Vincent Plagnol
}
