


qbetabinom <- function(p, size, phi, prob) {  ##parameterize with phi, prob instead of a,b
  a <- prob*(1-phi)/phi
  b <- (1-prob)*(1-phi)/phi
  return( qbetabinom.ab (p = p, size = size, shape1 = a, shape2 = b)  )
}


qbetabinom.ab <- function (p, size, shape1, shape2)  {
  my.p <- dbetabinom.ab(x = 0:size, , size = size, shape1 = shape1, shape2 = shape2)
  cs <- cumsum(my.p)
  above <- which(cs > p)[1]
  below <- above - 1
  
  if (below > 0) {
    return (below + (p - cs[below])/my.p[ above ])
  } else {
    return (p/my.p[1])
  }
}



viterbi.hmm <- function(transitions, loglikelihood, positions) {
  if ( nrow(transitions) != ncol(transitions) ) stop("Transition matrix is not square")
  nstates <- nrow(transitions)
  nobs <- nrow(loglikelihood)

  res <- .Call("C_hmm", nstates, nobs, transitions, loglikelihood, positions, PACKAGE = 'ExomeDepth')
  dimnames(res[[2]])[[2]] <- c('start.p', 'end.p', 'type', 'nexons')
  res[[2]] <- as.data.frame(res[[2]])
  names(res) <- c('Viterbi.path', 'calls')

  return(res)
}



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
      my.sim1 <- rbetabinom.ab( n = 2000, size = size, shape1 = my.alt.a, shape2 = my.alt.b)/size
                             
      my.BF <- dbeta(x = my.sim1, shape1 = my.alt.a, shape2 = my.alt.b, log = TRUE) -  dbeta(x = my.sim1, shape1 = my.a, shape2 = my.b, log = TRUE)   
      log10.my.BF <- log10(exp(1)) * my.BF
      my.res <- mean(log10.my.BF)
    }

    if (!limit) {
      my.pr <- dbetabinom.ab(x = 0:size, size = size, shape1 = my.alt.a, shape2 = my.alt.b, log = FALSE)
      my.BF <- dbetabinom.ab(x = 0:size, size = size, shape1 = my.alt.a, shape2 = my.alt.b, log = TRUE) - dbetabinom.ab(x = 0:size, size = size, shape1 = my.a, shape2 = my.b, log = TRUE)
      log10.my.BF <- log10(exp(1)) * my.BF
      my.res <- sum(my.pr * log10.my.BF)
    }
  }
  
  return(my.res)
}


