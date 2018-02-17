

#' Combine multiple samples to optimize the reference set in order to maximise
#' the power to detect CNV.
#'
#' The power to detect copy number variant (CNVs) from targeted sequence data
#' can be maximised if the most appropriate set of sequences is used as
#' reference. This function is designed to combine multiple reference exomes in
#' order to build the best reference set.
#'
#'
#' @param test.counts Read count data for the test sample (numeric, typically a
#' vector of integer values).
#' @param reference.counts Matrix of read count data for a set of additional
#' samples that can be used as a comparison point for the test sample.
#' @param bin.length Length (in bp) of each of the regions (often exons, but
#' not necessarily) that were used to compute the read count data (i.e. what is
#' provided in the argument test.counts of this function).  If not provided all
#' bins are assumed to have equal length.
#' @param n.bins.reduced This optimization function can be slow when applied
#' genome-wide. For the purpose of building the reference sample, it is not
#' necessary to use the full data. The number provided by this argument
#' specifies the number of regions (typically exons) that will be sub-sampled
#' (using a grid) to optimise the referenceset. I find that 10,000 is largely
#' sufficient for exome data.
#' @param data Defaults to NULL: A data frame of covariates that can be
#' included in the model.
#' @param formula Defaults to 'cbind(test, reference) ~ 1'. This formula will
#' be used to fit the read count data. Covariates present in the data frame
#' (for example GC content) can be included in the right hand side of the
#' equation'. If covariates are provided they must be provided as arguments (in
#' the data frame ``data'').
#' @param phi.bins Numeric integer (typically 1, 2, or 3) that specifies the
#' number of windows where the over-dispersion parameter phi can vary. It
#' defaults to 1, i.e. a single over-dispersion parameter, independently of
#' read depth.
#' @return \item{reference.choice }{character: list of samples selected as
#' optimum reference set.} \item{summary.stats}{A data frame summarizing the
#' output of this computation, including expected Bayes factor, Rs statistic
#' (see reference for explanation) for multiple choices of reference set.}


select.reference.set <- function(test.counts, reference.counts, bin.length = NULL, n.bins.reduced = 0, data = NULL, formula = 'cbind(test, reference) ~ 1', phi.bins = 1) {

  message('Optimization of the choice of aggregate reference set, this process can take some time')

  if (sum(test.counts > 2) < 5) {
    message('It looks like the test samples has only ', sum(test.counts > 2), ' bins with 2 or more reads. The coverage is too small to perform any meaningful inference so no likelihood will be computed.')
    my.res <- list(reference.choice = dimnames(reference.counts)[[2]][1], summary.stats = NULL)
    return(my.res)
  }
  

  
  if (class(reference.counts) != 'matrix') stop('The reference sequence count data must be provided as a matrix')
  if (nrow(reference.counts) != length(test.counts)) stop("The number of rows of the reference matrix must match the length of the test count data\n")
  if (is.null(bin.length)) bin.length <- rep(1, length(test.counts))
  
  
  n.ref.samples <- ncol(reference.counts)

  ############ select the subset of bins which will be used for the selection of the reference set
  total.counts <- apply(reference.counts, MARGIN = 1, FUN = sum) + test.counts
  my.quantiles <- quantile(total.counts [ which(total.counts > 30) ], prob = c(0.1, 0.9), na.rm = TRUE)
  
  selected <- which(total.counts > 30 &
                    bin.length >= as.numeric(quantile(bin.length, prob = 0.05, na.rm = TRUE)) & #I remove very small exons here, because they cause instability
                    bin.length <= as.numeric(quantile(bin.length, prob = 0.95, na.rm = TRUE)) &  # no large exons
                    total.counts < my.quantiles[2]) 
  if ( (n.bins.reduced < length(selected)) && (n.bins.reduced > 0) ) selected <- selected[ seq(1, length(selected), length(selected) / n.bins.reduced) ]


  test.counts <- test.counts[ selected ]
  #print(selected)



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

    my.mod <- new('ExomeDepth',
                  test = test.counts,
                  reference = reference,
                  formula = formula,
                  data = data,
                  phi.bins = phi.bins,
                  verbose = FALSE)
    
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
