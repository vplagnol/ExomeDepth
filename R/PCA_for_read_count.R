
correct.counts.using.PCA <- function( count.data, nPCs = 3, mask.exons = NULL ) {
  
  if (class(count.data) != 'matrix') {stop("The input to the PCA correction must be a matrix")}
  nsamples <- ncol(count.data)
  nexons <- nrow(count.data)

  message('Now applying the PCA, you provided a matrix with ', nexons, ' exons and ', nsamples, ' samples')
  
  norm.count <- count.data
  my.rsums <- rowSums(count.data)/1000
  for (i in 1:nsamples) norm.count[,i] <- norm.count[,i] / my.rsums[ i ]
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

