




select.reference.set <- function(test.counts, reference.counts, bin.length = NULL, n.bins.reduced = 0, data = NULL, formula = 'cbind(test, reference) ~ 1', phi.bins = 1) {

  message('Optimization of the choice of aggregate reference set')
  
  if (class(reference.counts) != 'matrix') stop('The reference sequence count data must be provided as a matrix')
  if (nrow(reference.counts) != length(test.counts)) stop("The number of rows of the reference matrix must match the length of the test count data\n")
  if (is.null(bin.length)) bin.length <- rep(1, length(test.counts))
  
  
  n.ref.samples <- ncol(reference.counts)

  ############ select the subset of bins which will be used for the selection of the reference set
  total.counts <- apply(reference.counts, MARGIN = 1, FUN = sum) + test.counts
  my.quantiles <- quantile(total.counts [ which(total.counts > 30) ], prob = c(0.1, 0.9))
  
  selected <- which(total.counts > 30 & bin.length > 0 & total.counts < my.quantiles[2])
  if ( (n.bins.reduced < length(selected)) && (n.bins.reduced > 0) ) selected <- selected[ seq(1, length(selected), length(selected) / n.bins.reduced) ]


  test.counts <- test.counts[ selected ]


  reference.counts <- reference.counts[ selected, , drop = FALSE ]
  bin.length <- bin.length[ selected]
  if (!is.null(data)) data <- data[ selected, ]
  n.bins <- length(selected)
  message('Number of selected bins: ', n.bins)

  ############### Now sort the data according to the correlation
  my.correlations <- apply(reference.counts, MARGIN = 2, FUN = function(x) {cor(x/(bin.length*sum(x)/10^6), test.counts/(bin.length*sum(test.counts)/10^6))})
  reference.counts <- reference.counts[, order(my.correlations, decreasing = TRUE), drop = FALSE]
  my.correlations <- my.correlations[ order(my.correlations, decreasing = TRUE) ]

  
  ########################################
  res.data.frame <- data.frame(ref.samples = dimnames(reference.counts)[[2]],
                               correlations = my.correlations,
                               expected.BF = NA,
                               phi = NA,
                               RatioSd = NA,
                               mean.p = NA,
                               median.depth = NA,
                               selected = FALSE)

  reference <- rep(0, n.bins)
  for (i in 1:n.ref.samples) {
    reference <- reference + reference.counts[,i]

    #print(head(data))
    #print(length(reference))
    #print(nrow(data))

    my.mod <- new('ExomeDepth',
                  test = test.counts,
                  reference = reference,
                  formula = formula,
                  data = data,
                  phi.bins = phi.bins)
    
    res.data.frame$phi[ i ] <- mean(my.mod@phi)
    res.data.frame$mean.p[ i ] <- mean(my.mod@expected)
    res.data.frame$median.depth[ i ] <- median(reference)
    res.data.frame$RatioSd[ i ] <-  mean(sqrt(1 + (test.counts + reference - 1)*my.mod@phi))

    if ( (i > 2) && (res.data.frame$mean.p[ i ] < 0.05)) break;
    
    ##########determine the expected proportion of reads assuming a deletion
    alt.odds <- res.data.frame$mean.p[ i ]/(1-res.data.frame$mean.p[ i ]) * 0.5
    alt.mean.p <- alt.odds/(1+alt.odds)

    res.data.frame$expected.BF[ i ] <- get.power.betabinom (size = round(res.data.frame$median.depth[ i ]),
                                                            my.phi = res.data.frame$phi[ i ],
                                                            my.p = res.data.frame$mean.p[ i ],
                                                            my.alt.p = alt.mean.p,
                                                            theory = FALSE)
  }

  my.max <- which.max( res.data.frame$expected.BF )
  res.data.frame$selected[ my.max ] <- TRUE
  my.res <- list(reference.choice = as.character(res.data.frame$ref.samples[ 1:my.max ]), summary.stats = res.data.frame)
                 
  return(my.res)
}
