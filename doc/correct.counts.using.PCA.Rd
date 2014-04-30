\name{correct.counts.using.PCA}
\alias{correct.counts.using.PCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Correct read count data using a PCA correction
}
\description{
  Applies a PCA correction to the matrix of read count to remove batch
  effects.
  Returns the residuals matrix after havnig applied a specified number
  of PC correction.
}
\usage{
correct.count.using.PCA(count.data, nPCs = 3)
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
}
\value{
  Returns a matrix with the same dimension as the input.
}
\author{
  Vincent Plagnol
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (count.data, nPCs = 3) 
{
    if (class(count.data) != "matrix") {
        stop("The input to the PCA correction must be a matrix")
    }
    nsamples <- ncol(count.data)
    nexons <- nrow(count.data)
    message("Now applying the PCA, you provided a matrix with", 
        nexons, " exons and", nsamples, "samples")
    norm.count <- count.data
    my.rsums <- rowSums(count.data)/1000
    for (i in 1:nsamples) norm.count[, i] <- norm.count[, i]/my.rsums[i]
    norm.count <- t(as.matrix(norm.count))
    centers <- colMeans(norm.count)
    for (i in 1:nexons) norm.count[, i] <- norm.count[, i] - 
        centers[i]
    my.pca <- prcomp(norm.count)
    PCA.mat <- my.pca$x[, 1:nPCs]
    reg.mat <- solve(t(PCA.mat) \%*\% PCA.mat) \%*\% t(PCA.mat)
    coeff.mat <- reg.mat \%*\% norm.count
    PCA.scores <- my.pca$x[, 1:nPCs] \%*\% coeff.mat
    norm.count <- norm.count - PCA.scores
    for (i in 1:nexons) norm.count[, i] <- round(pmax(0, my.rsums[i] * 
        (norm.count[, i] + centers[i])))
    norm.count <- t(norm.count)
    return(norm.count)
  }
}
