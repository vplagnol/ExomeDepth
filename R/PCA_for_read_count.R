norm.count <- count.data
my.rsums <- rowSums(count.data)/1000
for (i in 1:nsamples) norm.count[,i] <- norm.count[,i] / my.rsums[ i ]
norm.count <- t(as.matrix(norm.count))

############# Now prepare PCA
centers <-  colMeans(norm.count)
for (i in 1:nexons) norm.count[,i] <- norm.count[,i] - centers[ i ]
my.pca <- prcomp(norm.count)



pdf('fig/PCA.pdf')
plot (x = my.pca$x[,1], my.pca$x[,2], pch = '+')
plot (x = my.pca$x[,3], my.pca$x[,4], pch = '+')
dev.off()

PCA.mat <- my.pca$x[,1:nPCs]



reg.mat<-solve(t(PCA.mat) %*% PCA.mat)%*%t(PCA.mat)
coeff.mat<-reg.mat %*% norm.count

PCA.scores <- my.pca$x[, 1:nPCs] %*% coeff.mat

norm.count <- norm.count - PCA.scores
for (i in 1:nexons) norm.count[,i] <-  round(pmax(0, my.rsums[ i ]*(norm.count[,i] + centers[ i ])))
norm.count <- t(norm.count)
