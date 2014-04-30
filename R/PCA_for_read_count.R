
correct.counts.using.PCA <- function( count.data, nPCs = 3 ) {
  
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
  for (i in 1:nexons) norm.count[,i] <- norm.count[,i] - centers[ i ]
  my.pca <- prcomp(norm.count)
  
  PCA.mat <- my.pca$x[,1:nPCs]
  
  reg.mat<-solve(t(PCA.mat) %*% PCA.mat)%*%t(PCA.mat)
  coeff.mat<-reg.mat %*% norm.count
  PCA.scores <- my.pca$x[, 1:nPCs] %*% coeff.mat
  
  norm.count <- norm.count - PCA.scores
  for (i in 1:nexons) norm.count[,i] <-  round(pmax(0, my.rsums[ i ]*(norm.count[,i] + centers[ i ])))
  norm.count <- t(norm.count)
  return(norm.count)
}

