#### compute the power for a beta-binomial test


get.new.mean <- function(old.mean) {
  intercept <- log(old.mean/(1-old.mean))
  my.nu <- intercept + log(0.5)
  new.mean <- exp(my.nu)/(1+exp(my.nu))    
  return(new.mean)
}



get.power.betabin <- function (size, my.phi, my.p, my.alt.p, theory = FALSE, frequentist = FALSE, limit = FALSE) {
  require(VGAM)

  #my.alt.phi <- my.phi*my.alt.p/my.p
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
      my.sim1 <- rbetabin.ab( n = 2000, size = size, shape1 = my.alt.a, shape2 = my.alt.b)/size
                             
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



get.power.betabin.old <- function (size, my.phi, my.p, my.alt.p, theory = FALSE, frequentist = FALSE) {
  require(VGAM)
  
  if (!theory) {
    my.a <- my.p*(1-my.phi)/my.phi
    my.b <- (1-my.p)*(1-my.phi)/my.phi
    
    my.alt.a <- my.alt.p*(1-my.phi)/my.phi
    my.alt.b <- (1-my.alt.p)*(1-my.phi)/my.phi

    my.pr <- dbetabinom.ab(x = 0:size, size = size, shape1 = my.alt.a, shape2 = my.alt.b, log=FALSE)

    if (!frequentist) {
      my.BF <-  dbetabinom.ab(x = 0:size, size = size, shape1 = my.alt.a, shape2 = my.alt.b, log=TRUE) - dbetabinom.ab(x = 0:size, size = size, shape1 = my.a, shape2 = my.b, log=TRUE)
      log10.my.BF <- log10(exp(1))*my.BF
      return(sum(my.pr*log10.my.BF))
    }

    if (frequentist) {
      my.threshold <- qbetabin.ab(p = 10^-5, size = size, shape1 = my.a, shape = my.b)
      
      gap.should <- 10^-5 - pbetabin.ab(q = my.threshold-1, size = size, shape1 = my.a, shape2 = my.b)
      gap.is <- pbetabin.ab(q = my.threshold, size = size, shape1 = my.a, shape2 = my.b) - pbetabin.ab(q = my.threshold-1, size = size, shape1 = my.a, shape2 = my.b)
      my.power <-  pbetabin.ab(q = my.threshold-1, size = size, shape1 = my.alt.a, shape2 = my.alt.b) + (gap.should/gap.is)*dbetabinom.ab(x = my.threshold, size = size, shape1 = my.alt.a, shape2 = my.alt.b)
      return(my.power)
    }
    
  }
  

  ################
  if (theory) {

    if (!frequentist) {
      my.pr <- dbinom(x = 0:size, size = size, prob = my.alt.p, log=FALSE)
      my.BF <-   dbinom(x = 0:size, size = size, prob = my.alt.p, log=TRUE) - dbinom(x = 0:size, size = size, prob = my.p, log=TRUE)
      log10.my.BF <- log10(exp(1))*my.BF
      return(sum(my.pr*log10.my.BF))
    }
    
    if (frequentist) {
      my.threshold <- qbinom(p = 10^-5, size = size, prob = my.p)

      gap.should <- 10^-5 - pbinom(q = my.threshold-1, size = size, prob = my.p)
      gap.is <- pbinom(q = my.threshold, size = size, prob = my.p) - pbinom(q = my.threshold-1, size = size, prob = my.p)
      my.power <- pbinom(q = my.threshold-1, size = size, prob = my.alt.p) + (gap.should/gap.is)*dbinom(x = my.threshold, size = size, prob = my.alt.p)
      
      return (my.power)
    }
  }
  
  return(NA)
}







############### find the best reference set
exodepth.choice.reference <- function( my.data, my.test, excluded.samples, my.ref = NA, extra.computations = FALSE, all.rows = FALSE) {

  require(VGAM)
  require(aod)

  if (is.na(my.ref[1]))  my.ref <- grep (pattern = '*.bam', names(data), value = TRUE)   ##arg 3, may include the test sample
  
#####################  remove the exons that do not contribute
  my.sum.row <- apply(as.matrix(my.data[, unique(c(my.ref, my.test))]), MAR = 1, FUN = sum)
  data.loc <- my.data[ (my.sum.row > 10) & (my.data[, my.test] > 0) , ]

  
########## create the main data frame
  my.data.frame <- data.frame(samples = subset(my.ref, my.ref != my.test & ! (my.ref %in% excluded.samples) ),
                              mean.p = NA,
                              phi = NA,
                              phi.no.GC = NA,
                              expected.BF = NA,
                              expected.BF.no.GC = NA,
                              theory.BF = NA,
                              typical.depth = NA,
                              dispersion.parameter = NA,
                              cor.FPKM.single = NA,
                              cor.FPKM = NA,
                              choice = FALSE)

  data.loc$test <-  data.loc[, my.test]

  
################# compute basic correlations at the FPKM level
  for (i in 1:nrow(my.data.frame)) {
    my.s <- as.character(my.data.frame$samples[i])
    my.data.frame$cor.FPKM.single[i] <-  cor (data.loc$test/ (data.loc$length*sum(data.loc$test)/10^6), data.loc[, my.s] / (data.loc$length*sum(data.loc[, my.s])/10^6) )
  }
  my.data.frame <- my.data.frame[ order(my.data.frame$cor.FPKM.single, decreasing = TRUE), ]
  
##################

  
  data.loc$total.counts <- data.loc$test
  expected.BF.best <- -1
  for (i in 1:nrow(my.data.frame)) {
    my.s <- as.character(my.data.frame$samples[i])
    message(i, ' ', my.s, '\n')
    data.loc$total.counts <- data.loc$total.counts +  data.loc[, my.s]
    data.loc$failure <- data.loc$total.counts - data.loc$test

    expected.depth.ref <- mean(  data.loc$total.counts )


    ############### key fitting of parameters
    mod <- betabin( data = data.loc, formula = cbind(test, failure) ~ GC, random = ~ 1, link = 'logit')
                                        #eta <- mod@param[[ '(Intercept)' ]]
    
    eta <- mod@param[['GC']] * data.loc$GC + mod@param[['(Intercept)']]
    mean.p <- mean(exp(eta)/(1 + exp(eta) ))  ##average p, apply the inverse log-link      
    
    phi <-  mod@param[[ 'phi.(Intercept)']]
    ###########################

    my.n <- floor(median(data.loc$total.counts))    
    my.data.frame$mean.p[i] <- mean.p
    my.data.frame$phi[i] <- phi
    my.data.frame$typical.depth[i] <- my.n
    my.data.frame$dispersion.parameter[ i ] <- 1 + expected.depth.ref * phi
    my.data.frame$cor.FPKM[ i ] <- cor (data.loc$test/ (data.loc$length*sum(data.loc$test)/10^6), data.loc$failure / (data.loc$length*sum(data.loc$failure)/10^6) )
    
  ############ parameters under the alternate hypothesis
    alt.odds <- mean.p/(1-mean.p) * 0.5
    alt.mean.p <- alt.odds/(1+alt.odds)

    expected.BF <- get.power.betabin (size = my.n,
                                      my.phi = phi,
                                      my.p = mean.p,
                                      my.alt.p = alt.mean.p,
                                      theory = FALSE)
    
    my.data.frame$expected.BF[ i ] <-  expected.BF


    
    if (extra.computations) {

############## theoretical power
      my.data.frame$theory.BF[ i ] <- get.power.betabin (size = my.n,
                                                         my.phi = -1,
                                                         my.p = mean.p,
                                                         my.alt.p = alt.mean.p,
                                                         theory = TRUE)

############## Now without the GC content
      mod <- betabin( data = data.loc, formula = cbind(test, failure) ~ 1, random = ~ 1, link = 'logit')
      eta <- mod@param[['(Intercept)']]
      mean.p <- mean(exp(eta)/(1 + exp(eta) ))  ##average p, apply the inverse log-link
      alt.odds <- mean.p/(1-mean.p) * 0.5
      alt.mean.p <- alt.odds/(1+alt.odds)
      
      phi.no.GC <-  mod@param[[ 'phi.(Intercept)']]
      my.data.frame$phi.no.GC[i] <- phi.no.GC

      
      my.data.frame$expected.BF.no.GC[ i ] <- get.power.betabin (size = my.n,
                                                                 my.phi = phi.no.GC,
                                                                 my.p = mean.p,
                                                                 my.alt.p = alt.mean.p,
                                                                 theory = FALSE)
    }

###############
    message("Current read depth: ", signif(my.n, 3), ' and mean p: ', signif(mean.p, 3), ', expected log10BF: ', expected.BF)
    if ( (mean.p < 0.05) && (!all.rows) )  {print("Min p reached"); break;}
    
  }

  ##### returns the best guess
  my.max <- which.max(my.data.frame$expected.BF )
  my.data.frame$choice [ my.max ] <- TRUE
  res <- list(summary = my.data.frame, ref.samples = as.character(my.data.frame$samples[ 1:my.max ]))
  
}


