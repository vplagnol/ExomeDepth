#' Compute the residuals of the read count data after correcting for the first
#' principal components.
#'
#' Applies a PCA correction to the matrix of read count to remove batch
#' effects.  Returns the residuals matrix after having applied a correction for
#' a number of principal components specified by the user.
#'
#' Batch effects are unavoidable for most application of high throughput
#' sequencing technologies. ExomeDepth attempts to correct for this by building
#' an aggregate reference set that matches as closely as possible the test
#' sample. However, this may not be sufficient in all cases. This function
#' attempts to improve this situation by computing the sample covariance
#' matrix, computing the first principal components and taking the residuals of
#' the count data after having fitted these first few components.  This
#' function is relatively experimental so far, and is likely to work only if
#' the number of sample is relatively large ( > 20, typically, even though the
#' appropriate number is not clear). Performances on exome datasets have not
#' yet been tested thoroughly but it is best suited for targeted sequencing
#' experiments (i.e. smaller gene panels) for which the sample size is likely
#' to be larger and there should be no performance issue associcated with the
#' computation of the covariance matrix.  The \emph{mask.exons} argument allows
#' to mask exons for which one expects recurrent CNVs to be present. The
#' rationale being that if these deletions are common, they may contribute
#' largely to the first few principal components and remove the signal. Common
#' CNVs should also probably be removed because one typically attempts to
#' remove technical artifacts rather than actual variability in copy number.
#'
#' @param count.data A matrix of read count, with samples as columns and exons
#' as rows.
#' @param nPCs Number of PCs that must be applied when computing the residuals
#' of the read count data.
#' @param mask.exons A logical vector that indicates which exons should be
#' masked when the sample covariance matrix is computed (see details for more
#' information).
#' @return Returns a matrix of corrected read count data of the same dimension
#' as the input.
#' @author Vincent Plagnol



correct.counts.using.PCA <- function( count.data, nPCs = 3, mask.exons = NULL ) {
  
  if (class(count.data) != 'matrix') {stop("The input to the PCA correction must be a matrix")}
  nsamples <- ncol(count.data)
  nexons <- nrow(count.data)

  message('Now applying the PCA, you provided a matrix with ', nexons, ' exons and ', nsamples, ' samples')
  
  norm.count <- count.data
  my.rsums <- rowSums(count.data)/1000  ##this is used to normalize for the read depth, and I normalize it by 1,000 reads
  for (i in 1:nsamples) norm.count[,i] <- norm.count[,i] / max(1, my.rsums[ i ])
  norm.count <- t(as.matrix(norm.count))
  
############# Now prepare PCA
  centers <-  colMeans(norm.count)
  good.depth <- apply(MAR = 2, norm.count, FUN = sd) > 2  ##here I want to remove exons that do not have any variability (all 0s usually)
  
  for (i in 1:nexons) norm.count[,i] <- norm.count[,i] - centers[ i ]
  if (!is.null(mask.exons)) {
    if (class(mask.exons) != 'logical') stop('The mask exons argument must be a logical vector')
    if (length(mask.exons) != nexons) stop('The length of the mask exons argument does not match the number of exons')
    message('You are masking ', sum(mask.exons), ' exons in the PCA computation step')
    my.pca <- prcomp(norm.count[, !mask.exons & good.depth])
  } else {
    my.pca <- prcomp(norm.count[, good.depth ])
  }
  
  PCA.mat <- my.pca$x[,1:nPCs]
  
  reg.mat<-solve(t(PCA.mat) %*% PCA.mat)%*%t(PCA.mat)
  coeff.mat<-reg.mat %*% norm.count
  PCA.scores <- my.pca$x[, 1:nPCs] %*% coeff.mat
  
  norm.count <- norm.count - PCA.scores
  for (i in 1:nexons) norm.count[,i] <-  round(pmax(0, my.rsums[ i ]*(norm.count[,i] + centers[ i ])))
  norm.count <- t(norm.count)
  return(norm.count)
}

