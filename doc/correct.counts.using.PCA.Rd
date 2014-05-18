\name{correct.counts.using.PCA}
\alias{correct.counts.using.PCA}
\title{
  Compute the residuals of the read count data after correcting for the
  first principal components.
}
\description{
  Applies a PCA correction to the matrix of read count to remove batch
  effects.
  Returns the residuals matrix after having applied a correction for a
  number of principal components specified by the user.
}
\usage{
correct.count.using.PCA(count.data, nPCs = 3, mask.exons = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count.data}{
    A matrix of read count, with samples as columns and exons as rows.
}
  \item{nPCs}{
    Number of PCs that must be applied when computing the residuals of
    the read count data.
  }
  \item{mask.exons}{
  A logical vector that indicates which exons should be masked when the
  sample covariance matrix is computed (see details for more information)}
}
\value{
  Returns a matrix of corrected read count data of the same dimension as the input.
}
\details{
  Batch effects are unavoidable for most application of high throughput
  sequencing technologies. ExomeDepth attempts to correct for this by
  building an aggregate reference set that matches as closely as
  possible the test sample. However, this may not be sufficient in all
  cases. This function attempts to improve this situation by computing
  the sample covariance matrix, computing the first principal components
  and taking the residuals of the count data after having fitted these
  first few components.
  This function is relatively experimental so far, and is only likely to
  work if the number of sample is relatively large ( > 20, typically,
  even though the appropriate number is not clear). Performances on
  exome datasets have not yet been tested thoroughly but it is best
  suited for targeted sequencing experiments (i.e. smaller gene panels)
  for which the sample size is likely to be larger and there should be
  no performance issue associcated with the computation of the
  covariance matrix.
  The \texttt{mask.exons} argument allows to mask exons for which one
  expects recurrent CNVs to be present. The rationale being that if these
  deletions are common, they may contribute largely to the first few
  principal components and remove the signal. Common CNVs should also
  probably be removed because one typically attempts to remove technical
  artifacts rather than actual variability in copy number.
}
\author{
  Vincent Plagnol
}
