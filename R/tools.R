#' Quantile for betabin function
#'
#' Quantile function for the betabinomal distribution using the p/phi
#' parameterisation.
#'
#' Filling a gap in the VGAM package.
#'
#' @param p Point of the distribution from which one is looking for the
#' quantile
#' @param size Sample size of the random variable
#' @param phi Over-dispersion parameter
#' @param prob Mean probability of the binomial distribution
#' @return A real number corresponding to the quantile p.
#' @author Vincent Plagnol
#' @seealso VGAM R package.

qbetabinom <- function(p, size, phi, prob) {  ##parameterize with phi, prob instead of a,b
  a <- prob*(1-phi)/phi
  b <- (1-prob)*(1-phi)/phi
  return( qbetabinom.ab (p = p, size = size, shape1 = a, shape2 = b)  )
}



#' Quantile function for the beta-binomial distribution
#'
#' Standard qbetabinomial.ab function which is missing from the VGAM package.
#'
#'
#' @param p Mean value of the beta-binomial distribution.
#' @param size Size of the beta-binomial.
#' @param shape1 First parameter of the beta distribution for p.
#' @param shape2 Second parameter of the beta distribution for p.
#' @return A quantile of the distribution.
#' @seealso \code{VGAM} package.

qbetabinom.ab <- function (p, size, shape1, shape2)  {
  my.p <- VGAM::dbetabinom.ab(x = 0:size, , size = size, shape1 = shape1, shape2 = shape2)
  cs <- cumsum(my.p)
  above <- which(cs > p)[1]
  below <- above - 1

  if (below > 0) {
    return (below + (p - cs[below])/my.p[ above ])
  } else {
    return (p/my.p[1])
  }
}


#' Computes the Viterbi path for a hidden markov model
#'
#' Estimates the most likely path for a hidden Markov Chain using the maximum
#' likelihood Viterbi algorithm.
#'
#' Standard forward-backward Viterbi algorithm using a precomputed matrix of
#' likelihoods.
#'
#' @param transitions Transition matrix
#' @param loglikelihood numeric matrix containing the loglikelihood of the data
#' under the possible states
#' @param positions Positions of the exons
#' @param expected.CNV.length Expected length of CNV calls, which has an impact
#' on the transition matrix between CNV states.
#' @return \item{comp1 }{Description of 'comp1'} \item{comp2 }{Description of
#' 'comp2'}


viterbi.hmm <- function(transitions, loglikelihood, positions, expected.CNV.length) {
  if ( nrow(transitions) != ncol(transitions) ) stop("Transition matrix is not square")
  if ( length(positions) != nrow(loglikelihood) ) {
    stop("The number of positions are not matching the number of rows of the likelihood matrix", length(positions), " and ", nrow(loglikelihood))
  }

  nstates <- nrow(transitions)
  nobs <- nrow(loglikelihood)

  res <- .Call("C_hmm", nstates, nobs, transitions, loglikelihood, positions, as.double(expected.CNV.length), PACKAGE = 'ExomeDepth')
  dimnames(res[[2]])[[2]] <- c('start.p', 'end.p', 'type', 'nexons')
  res[[2]] <- as.data.frame(res[[2]])
  names(res) <- c('Viterbi.path', 'calls')

  return(res)
}



#' Estimate the power to compare two beta-binomial distributions.
#'
#' A power study useful in the context of ExomeDepth.
#'
#'
#' @param size Number of samples from the beta-binomial distribution.
#' @param my.phi Over-dispersion parameter.
#' @param my.p Expected p under the null.
#' @param my.alt.p Expected p under the alternative.
#' @param theory \code{logical}, should a theoretical limit (large sample size)
#' be used? Defaults to FALSE.
#' @param frequentist \code{logical}, should a frequentist version be used?
#' Defaults to FALSE.
#' @param limit \code{logical}, should another large sample size limit be used?
#' Defaults to FALSE.
#' @return An expected Bayes factor.


get.power.betabinom <- function (size, my.phi, my.p, my.alt.p, theory = FALSE, frequentist = FALSE, limit = FALSE) {

  my.alt.phi <- my.phi

  my.a <- my.p * (1 - my.phi)/my.phi
  my.b <- (1 - my.p) * (1 - my.phi)/my.phi
  my.alt.a <- my.alt.p * (1 - my.alt.phi)/my.alt.phi
  my.alt.b <- (1 - my.alt.p) * (1 - my.alt.phi)/my.alt.phi

  if (theory) {  ## binomial case
    my.pr <- dbinom(x = 0:size, size = size, prob = my.alt.p, log=FALSE)
    my.BF <-   dbinom(x = 0:size, size = size, prob = my.alt.p, log=TRUE) - dbinom(x = 0:size, size = size, prob = my.p, log=TRUE)
    log10.my.BF <- log10(exp(1))*my.BF
    my.res <- sum(my.pr * log10.my.BF)
  }


  if (!theory) {  ##betabinomial case



    if (limit) { ##### beta approx

      #my.sim1 <-  rbeta(n = 2000,  shape1 = my.alt.a, shape2 = my.alt.b)
      my.sim1 <- VGAM::rbetabinom.ab( n = 2000, size = size, shape1 = my.alt.a, shape2 = my.alt.b)/size

      my.BF <- dbeta(x = my.sim1, shape1 = my.alt.a, shape2 = my.alt.b, log = TRUE) -  dbeta(x = my.sim1, shape1 = my.a, shape2 = my.b, log = TRUE)
      log10.my.BF <- log10(exp(1)) * my.BF
      my.res <- mean(log10.my.BF)
    }

    if (!limit) {
      my.pr <- VGAM::dbetabinom.ab(x = 0:size, size = size, shape1 = my.alt.a, shape2 = my.alt.b, log = FALSE)
      my.BF <- VGAM::dbetabinom.ab(x = 0:size, size = size, shape1 = my.alt.a, shape2 = my.alt.b, log = TRUE) - VGAM::dbetabinom.ab(x = 0:size, size = size, shape1 = my.a, shape2 = my.b, log = TRUE)
      log10.my.BF <- log10(exp(1)) * my.BF
      my.res <- sum(my.pr * log10.my.BF)
    }
  }

  return(my.res)
}


