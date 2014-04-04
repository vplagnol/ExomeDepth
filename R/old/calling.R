
######################## estimate the likelihood for different copy numbers
get.likelihood.exodepth <- function(ndata, formula = 'cbind(test, reference) ~ 1') {
  require(aod)
  
  for (my.na in c('test', 'reference')) {
    if (! my.na %in% names(ndata)) stop("Missing column ", my.na, " in the ndata data.frame")
  }
  
  mod <- betabin( data = ndata, formula = as.formula(formula), random = ~ 1, link = 'logit', hessian = FALSE)
  ndata$expected <- aod::fitted(mod)
  phi <-  mod@param[[ 'phi.(Intercept)']]
  
################################# Now get the matrix of likelihood
  my.likelihood <- .Call("get_loglike_matrix", phi = phi, expected = ndata$expected, total = as.integer(ndata$test + ndata$reference), observed = as.integer(ndata$test), PACKAGE = 'ExoDepth')
  ndata$loglike.del <- my.likelihood[,1]
  ndata$loglike <- my.likelihood[,2]
  ndata$loglike.dup <- my.likelihood[,3]
  ndata$call <- ifelse ( ndata$loglike.del >  ndata$loglike + 3, 'del',  ifelse ( ndata$loglike.dup >  ndata$loglike + 3, 'dup', NA))
  
  return(ndata)
}



##############################################
##only test and reference are required

call.CNVs <- function (test, reference, prior = 10^(-4), GC, chromosome, start, end, single.exon.output.file) {
  require(aod)
  
  if (length(test) != length(reference)) stop("Length of test and reference do not match")

  data <- data.frame (test = test, reference = reference, total = test + reference)
  
  if (!missing(GC)) {
    if ( sum(is.na(GC)) > 0 || (max(GC) > 1) || min(GC) < 0 ) stop("If provided, GC must be between 0 and 1 and have no missing data")
    if (length(GC) != length(test)) stop("Length of GC content do not match")
    data$GC <- GC
    my.formula <- 'cbind(test, reference) ~ GC'
    message("GC content incorporated into the model")
  } else {
    my.formula <- 'cbind(test, reference) ~ 1'
  }

  if (!missing(chromosome)) {
    if (length(chromosome) != length(test))  stop("Length of chromosome vector does not match")
    data$chromosome <- chromosome
  }
  
  if (!missing(start)) {
    if (length(start) != length(test))  stop("Length of chromosome vector does not match")
    data$start <- start
  }
  if (!missing(end)) {
    data$end <- end
    if (length(end) != length(test))  stop("Length of chromosome vector does not match")
  }


  data <- get.likelihood.exodepth (data, formula = my.formula)   ######## estimate the key parameters

  if (!missing(single.exon.output.file)) {
    interesting <- subset(data, !is.na(data$call))
    key.cols <- c('chromosome', 'start', 'end')
    interesting <- interesting[, c(key.cols, subset(names(interesting), ! names(interesting) %in% key.cols)) ]
                               
    write.table(x = format(interesting, digits = 3) ,
                file = single.exon.output.file,
                row.names = FALSE,
                quote = FALSE,
                sep = '\t')
  }

  
  if (!missing(chromosome)) {
    my.breaks <- which(diff(as.numeric(data$chromosome)) != 0) + 1
    data$loglike.del[ my.breaks ] <- - Inf
    data$loglike.dup[ my.breaks ] <- - Inf
  }
  


  proba <- as.matrix(data[, c('loglike', 'loglike.del', 'loglike.dup')])
  transitions <- matrix(nrow = 3, ncol = 3,
                        c( 1. - prior, prior/2., prior/2.,
                          0.5, 0.5, 0.,
                          0.5, 0, 0.5),
                        byrow = TRUE)
  positions <- (data$start + data$end)/2
  my.calls <- viterbi.hmm (transitions, loglikelihood = proba, positions = positions)

  if (nrow(my.calls$calls) > 0) {
################ outputs the data
    if (!missing(start)) my.calls$calls$start <- data$start[ my.calls$calls$start.p ] else my.calls$calls$start <- my.calls$calls$start.p
    if (!missing(end)) my.calls$calls$end <- data$end[ my.calls$calls$end.p ] else my.calls$calls$end <-  my.calls$calls$end.p 
    if (!missing(chromosome)) my.calls$calls$chromosome <- data$chromosome[ my.calls$calls$start.p ] else  my.calls$calls$chromosome <- ''
  
    my.calls$calls$id <- paste('chr', my.calls$calls$chromosome, ':',  my.calls$calls$start, '-',  my.calls$calls$end, sep = '')
    my.calls$calls$type <- c('deletion', 'duplication')[ my.calls$calls$type ]
  
########## make things pretty
    my.calls$calls$BF <- NA
    my.calls$calls$reads.expected <- NA
    my.calls$calls$reads.observed <- NA
  
    for (ir in 1:nrow(my.calls$calls)) {
      
      if (my.calls$calls$type[ir] == 'duplication') my.calls$calls$BF[ir] <-  sum(data$loglike.dup [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] - data$loglike [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ])
    
  
      if (my.calls$calls$type[ir] == 'deletion') my.calls$calls$BF[ir] <-  sum(data$loglike.del [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] - data$loglike [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ])
      
      my.calls$calls$reads.expected[ ir ] <-  sum( data$total [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] * data$expected [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ ir ] ])
      my.calls$calls$reads.observed[ ir ] <-  sum( data$test [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] )
      
    }
    
    my.calls$calls$reads.expected <- as.integer( my.calls$calls$reads.expected)
    my.calls$calls$reads.ratio <-  signif(my.calls$calls$reads.observed / my.calls$calls$reads.expected, 3)
    my.calls$calls$BF <- signif( log10(exp(1))*my.calls$calls$BF, 3)
  }
    
  return(my.calls)
}






