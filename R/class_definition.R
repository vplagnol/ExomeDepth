
#' Class `ExomeDepth`
#'
#' A class to hold the read count data that is used by ExomeDepth to call CNVs.
#'
#'
#' @name ExomeDepth-class
#' @md
#' @aliases ExomeDepth-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' `new("ExomeDepth", data = NULL, test, reference, formula = 'cbind(test,
#' reference) ~ 1', subset.for.speed = NULL)`.  `data` is optional and is
#' only used if the `formula` argument refers to covariates (in which case
#' these covariates must be included in the data frame).  `test` and
#' `reference` refer to the read count data for the test and reference
#' samples.
#'
#' Critically, it is not required to store the positions of the DNA fragments
#' that led to the test and reference counts. That is only required for the
#' function `TestCNV`. If this is of use, a GRanges object can be provided
#' using the argument `positions`.
#'
#' Creating a ExomeDepth object will automatically fit the
#' beta-binomial model (using routines from the `aod` package) and compute
#' the likelihood for the three copy number states (normal, deletion and
#' duplication).
#'
#' @seealso `?select.reference.set` `?CallCNVs`
#' @references A robust model for read count data in exome sequencing experiments and implications for copy number variant calling, Plagnol et al 2012
#' @keywords classes
#' @examples
#' showClass("ExomeDepth")


setClass("ExomeDepth",
         representation(test = "numeric",
                        reference = "numeric",
                        formula = "character",
                        expected = "numeric",
                        phi = "numeric",
                        likelihood = "matrix",
                        annotations = "data.frame",
                        positions = "GRanges",
                        CNV.calls = "data.frame"))




#' @title ExomeDepth initialization tool
#' @description Builds an exomeDepth object from test and reference vectors
#' @param .Object ExomeDepth object
#' @param data Data frame containing potential covariates.
#' @param test Numeric, vector of counts for the test sample.
#' @param reference Numeric, vector of counts for the reference sample.
#' @param formula Linear model to be used when fitting the data.
#' @param phi.bins Numeric, defaults to 1. Number of discrete bins for the over-dispersion parameter phi, depending on read depth.
#' Do not modify this parameter for the standard use of ExomeDepth.
#' @param prop.tumor Numeric, defaults to 1. For the somatic variant calling, this assesses the proportion of the test sample data originating from the tumour.
#' Do not modify this parameter for the standard use of ExomeDepth.
#' @param subset.for.speed Numeric, defaults to NULL. If non-null, this sets the number of data points to be used for an accelerated fit of the data.
#' @param positions Optional GRanges argument specifying the positions of the exons (or DNA regions) where the reads were counted for test and reference.
#' @param verbose Logical, controls the output level.
#' @return An ExomeDepth object, which contains the CNV calls after running a Viterbi algorithm.
#' @examples
#'
#' data(ExomeCount)  #pick an example count file
#' small_count <- ExomeCount[1:100, ]  #reduce the size for speedy computations
#'
#' ## remove exons without data below
#' small_count <- small_count[ small_count$Exome2 + small_count$Exome3 > 0, ]
#'
#' example_object <- new('ExomeDepth', test = small_count$Exome2,
#'                                     reference = small_count$Exome3,
#'                                     formula = 'cbind(test, reference) ~ 1')
#' print(example_object)
#' print( mean(example_object@expected)) ## proportion of reads expected to match the test set



setMethod("initialize", "ExomeDepth", function(.Object,
                                               data = NULL,
                                               test,
                                               reference,
                                               formula = 'cbind(test, reference) ~ 1',
                                               phi.bins = 1,
                                               prop.tumor = 1,
                                               subset.for.speed = NULL,
                                               positions = GenomicRanges::GRanges(),
                                               verbose = TRUE) {
  if (length(test) != length(reference)) stop("Length of test and numeric must match")

  if (sum(test > 5) < 5) {
    message('It looks like the test samples has only ', sum(test > 5), ' bins with more than 5 reads. The coverage is too small to perform any meaningful inference so no likelihood will be computed.')
    return(.Object)
  }


  n.data.points <- length(test)
  if (is.null(data)) data <- data.frame(intercept = rep(1, length(test)))

  ### the stuff below will be used for the estimation of phi
  data$test <- test
  data$reference <- reference

  if (!is.null(subset.for.speed)) {

    if ( is.numeric(subset.for.speed) && (length(subset.for.speed) == 1)) {subset.for.speed <- seq(from = 1, to = nrow(data), by = floor( nrow(data) / subset.for.speed ) )}  ###deals with unexpected use of subset.for.speed
    subset.for.speed <- subset.for.speed[ subset.for.speed %in% 1:nrow(data) ] ##make sure we do not select non-existing rows
    data.for.fit <- data[ subset.for.speed, , drop = FALSE]
  } else {
    data.for.fit <- data
  }

  if (verbose) message('Now fitting the beta-binomial model on a data frame with ', nrow(data.for.fit), ' rows : this step can take a few minutes.')
    if (phi.bins == 1) {
      mod <- aod::betabin(data = data.for.fit, formula = as.formula(formula), random = ~ 1, link = 'logit', warnings = FALSE)
      .Object@phi <- rep(mod@param[[ 'phi.(Intercept)']], n.data.points)
    } else {
      if (!is.null(subset.for.speed)) {stop('Subset for speed option is not compatible with variable phi. This will be fixed later on but for now please adapt your code.')}
      ceiling.bin <- quantile(reference, probs = c( 0.85, 1) )
      bottom.bins <- seq(from = 0, to =  ceiling.bin[1], by =  ceiling.bin[1]/(phi.bins-1))
      complete.bins <- as.numeric(c(bottom.bins, ceiling.bin[2] + 1))
      data$depth.quant <- factor(sapply(reference, FUN = function(x) {sum (x >= complete.bins)}))

      ## a check
      my.tab <-  table(data$depth.quant)
      if (length(my.tab) != phi.bins) {
        stop('Binning did not happen properly')
      }

      mod <- aod::betabin (data = data.for.fit,
                           formula = as.formula(formula),
                           random = as.formula('~ depth.quant'),
                           link = 'logit',
                           warnings = FALSE)

      phi.estimates <-  as.numeric(mod@random.param)
      data$phi <-  phi.estimates[  data$depth.quant ]

      ## Now the linear interpolation for phi, the over-dispersion parameter
      fc <- approxfun (x = c(complete.bins[ 1:phi.bins] + complete.bins[ 2:(phi.bins+1)])/2, y =  phi.estimates, yleft = phi.estimates[1], yright = phi.estimates[phi.bins])
      data$phi.linear <- fc (reference)

      .Object@phi <- data$phi.linear
    }



  .Object@positions <- positions
  .Object@formula <- formula
  .Object@test <- test
  .Object@reference <- reference
  my.coeffs <- mod@fixed.param


  ## if positions were provided, check that they match the expected length
  if (length(.Object@positions) > 0) {
      if (length(.Object@positions) != length(test)) {
          stop("The provided genomic positions (GRanges object) and test count vector are not matching in length: ", length(.Object@positions), ", ", length(test))
      }
  }

  if (is.null(subset.for.speed)) {
      ## looks like a bug in model.matrix(mt, mfb, contrasts) in aod::fitted
      ## so suppressing warnings below
      .Object@expected <- suppressWarnings(aod::fitted(mod))
  } else {
    intercept <- my.coeffs[[ '(Intercept)' ]]
    .Object@expected <- rep(intercept, times = nrow(data))
    if (length(my.coeffs) > 1) {
      for (na in names(my.coeffs)[ -1 ]) {
        .Object@expected = .Object@expected + my.coeffs [[ na ]]*data[, na]
      }
    }
    .Object@expected <- exp(.Object@expected)/ (1 + exp(.Object@expected))
  }

  .Object@annotations <- data.frame()

  if (verbose) message('Now computing the likelihood for the different copy number states')
  if (prop.tumor < 1) message('Proportion of tumor DNA is ', prop.tumor)
  .Object@likelihood <- .Call("get_loglike_matrix",
                              phi = .Object@phi,
                              expected = .Object@expected,
                              total = as.integer(.Object@reference + .Object@test),
                              observed = as.integer(.Object@test),
                              mixture = prop.tumor)
  .Object
})


#############################################################################
if (!isGeneric("show")) {
  if (is.function("show"))
    fun <- show
  else fun <- function(object) standardGeneric("show")
  setGeneric("show", fun)
}

setMethod("show", "ExomeDepth", function(object) {
  cat('Number of data points: ', length(object@test), '\n')
  cat('Formula: ', object@formula, '\n')
  cat('Phi parameter (range if multiple values have been set): ', range(object@phi), '\n')
  if (ncol(object@likelihood) == 3) cat("Likelihood computed\n") else cat("Likelihood not computed\n")
})

#############################################################################


setGeneric("TestCNV", def = function(x, chromosome, start, end, type) standardGeneric('TestCNV'))

#' @title TestCNV
#' @md
#' @aliases TestCNV
#' @description Computes the Bayes Factor in favour of a CNV defined by position and type.
#' @param x ExomeDepth object
#' @param chromosome Character, chromosome name.
#' @param start Numeric, start of the tested CNV
#' @param end Numeric, end of the tested CNV
#' @param type Character, must be either `deletion` or `duplication`.
#' @return A single numeric value that is the log likelihood ratio in favour of a call present at this location.
#' @examples
#' data(ExomeCount)
#' ExomeCount <- ExomeCount[1:500,] ## small for the purpose of this test
#' ref_counts <- ExomeCount$Exome2 + ExomeCount$Exome3 + ExomeCount$Exome4
#'
#' ## creates a simple ExomeDepth object
#' ## Note that I include the positions here (GRanges format)
#' ## This is necessary for TestCNV to work
#' test_object <- new('ExomeDepth', test = ExomeCount$Exome1,
#'                                  reference = ref_counts,
#'                                  positions = ExomeCount)
#'
#' ## A positive control first
#' TestCNV (test_object, chromosome = 'chr1', start = 1387428, end = 1405539, type = 'deletion')
#'
#' ## And a region without evidence of call
#' TestCNV (test_object, chromosome = 'chr1', start = 987428, end = 1005539, type = 'deletion')
#'

setMethod("TestCNV", "ExomeDepth", function(x, chromosome, start, end, type) {
  if (! type %in% c('deletion', 'duplication')) stop("type must be either duplication or deletion\n")
  if (length(chromosome) != 1 || length(start) != 1 || length(end) != 1 || length(type) != 1) stop("The arguments chromosome, start, end and type must all be of length 1")
  if (is.factor(chromosome)) chromosome <- as.character(chromosome)
  if (!is.character(chromosome)) stop('The input chromosome must be a character or a factor')

  if (length(x@positions) == 0) stop("This function cannot be used if the position of the exons/DNA segments was not included in the ExomeDepth object")

  which.exons <- which(as.logical( (GenomicRanges::seqnames(x@positions) == chromosome) & (GenomicRanges::start(x@positions) >= start) & (GenomicRanges::end(x@positions) <= end)))

  if (type == 'deletion') log.ratio <- sum(x@likelihood[ which.exons, 1] - x@likelihood[ which.exons, 2])
  if (type == 'duplication') log.ratio <- sum(x@likelihood[ which.exons, 3] - x@likelihood[ which.exons, 2])
  return  (log.ratio)
})



setGeneric("CallCNVs",
           def = function(x, chromosome, start, end, name, transition.probability = 0.0001, expected.CNV.length = 50000) standardGeneric('CallCNVs'))


#' @title CallCNVs
#' @aliases CallCNVs
#' @description Call CNV data from an ExomeDepth object.
#' @md
#'
#' @details Must be called on an ExomeDepth object and fits a hidden Markov model to the read depth data with three hidden states (normal, deletion, duplication).
#' Likelihood data must have been pre-computed which should have been done by default when the ExomeDepth object was created.
#'

#' @param x An `ExomeDepth` object
#' @param chromosome Chromosome information for each exon (factor).
#' @param start Start (physical position) of each exon (numeric, must have the
#' same length as the chromosome argument).
#' @param end End (physical position) of each exome (numeric, must have the
#' same length as the chromosome argument).
#' @param name Name of each exon (character or factor).
#' @param transition.probability Transition probability of the hidden Markov
#' Chain from the normal copy number state to either a deletion or a
#' duplication. The default (0.0001) expect approximately 20 CNVs genome-wide.
#' @param expected.CNV.length The expectation for the length of a CNV. This
#' value factors into the Viterbi algorithm that is used to compte the
#' transition from one state to the next, which depends on the distance between
#' exons.
#' @return The same ExomeDepth object provided as input but with the slot
#' CNVcalls containing a data frame with the output of the calling.
#' @examples
#'
#' data(ExomeCount)
#' ExomeCount <- ExomeCount[1:500,] ## small for the purpose of this test
#' ref_counts <- ExomeCount$Exome2 + ExomeCount$Exome3 + ExomeCount$Exome4
#'
#' ## creates a simple ExomeDepth object
#' test_object <- new('ExomeDepth', test = ExomeCount$Exome1, reference = ref_counts)
#'
#' ## Call CNVs
#' called_object <- CallCNVs(x = test_object, transition.probability = 10^-4,
#'                          chromosome = GenomicRanges::seqnames(ExomeCount),
#'                          start = GenomicRanges::start(ExomeCount),
#'                          end = GenomicRanges::end(ExomeCount),
#'                          name = ExomeCount$names)
#'
#'
#' print(called_object@CNV.calls)
#'



setMethod("CallCNVs", "ExomeDepth", function( x, chromosome, start, end, name, transition.probability, expected.CNV.length) {

  if (length(x@phi) == 0) {
    message('The vector phi does not seem initialized. This may be because the read count is too low and the test vector cannot be processed. No calling will happen')
    x@CNV.calls <- data.frame()
    return(x)
  }

  if ( length(start) != length(chromosome) || length(end) != length(chromosome) || length(name) != length(chromosome) ) stop('Chromosome, start and end vector must have the same lengths.\n')
  if (nrow(x@likelihood) != length(chromosome) ) stop('The annotation vectors must have the same length as the data in the ExomeDepth x')

  ### Try to get the chromosome order right
  chr.names.used <- unique(as.character(chromosome))
  chr.levels <- c(as.character(seq(1, 22)),  chr.names.used[! chr.names.used %in% as.character(seq(1, 22)) ] )
  chr.levels <- chr.levels[ chr.levels %in% chr.names.used ]

  x@annotations <- data.frame(name = name, chromosome = factor(chromosome, levels = chr.levels), start = start, end = end)
  my.new.order <-  order(x@annotations$chromosome, 0.5*(x@annotations$start + x@annotations$end) )

  if (sum( my.new.order != 1:nrow(x@annotations) ) > 0) {
    message('Positions of exons seem non ordered, so ExomeDepth will reorder the data according to chromosome and position')
    x@test <- x@test[ my.new.order ]
    x@reference <- x@reference[ my.new.order ]
    x@annotations <- x@annotations[ my.new.order, ]
    x@likelihood <- x@likelihood[ my.new.order, ]
  }

  cor.test.reference <- cor(x@test, x@reference)
  message('Correlation between reference and tests count is ', signif(cor.test.reference, 5))
  message('To get meaningful result, this correlation should really be above 0.97. If this is not the case, consider the output of ExomeDepth as less reliable (i.e. most likely a high false positive rate)')

  total <- x@test + x@reference
  transitions <- matrix(nrow = 3, ncol = 3,
                        c( 1. - transition.probability, transition.probability/2., transition.probability/2.,
                          0.5, 0.5, 0.,
                          0.5, 0, 0.5),
                        byrow = TRUE)

  my.chromosomes <- unique(x@annotations$chromosome)

  final <- data.frame()

  shift <- 0
  for (chrom in my.chromosomes)  {  ##run Viterbi separately for each chromosome
    good.pos <- which (x@annotations$chromosome == chrom)
    loc.annotations <- x@annotations[  good.pos , ]
    loc.expected <- x@expected[ good.pos ]
    loc.test <- x@test[ good.pos ]
    loc.total <- total[ good.pos ]
    positions <- loc.annotations$start
    end.positions <- loc.annotations$end       #end position of targeted exons to be used when adding a dummy exon at the end of the chromosome.

    ##loc.likelihood <-  rbind(c(- Inf, 0, -Inf), x@likelihood[good.pos, c(2, 1, 3)]) ##add a dummy exon so that we start at cn = 2 (normal)
    loc.likelihood <- rbind(c(- Inf, 0, -Inf), x@likelihood[good.pos, c(2, 1, 3)],c(-100,0,-100)) ##update from Anna Fowler add a dummy exon so that we start at cn = 2 (normal) and a dummy exon at the end of the chromosome as well so that it ends at cn=2; for some reason I had to use -100 instead of -Inf

    my.calls <- viterbi.hmm (transitions,
                             loglikelihood = loc.likelihood,
                             positions = as.integer(c(positions[1] - 2*expected.CNV.length, positions,end.positions[length(end.positions)]+2*expected.CNV.length)),   #include position of new dummy exon
                             expected.CNV.length = expected.CNV.length)

    my.calls$calls$start.p <- my.calls$calls$start.p -1  ##remove the dummy exon, which has now served its purpose
    my.calls$calls$end.p <- my.calls$calls$end.p -1  ##remove the dummy exon, which has now served its purpose
    #loc.likelihood <- loc.likelihood[ -1, c(2,1, 3) ]  ##remove the dummy exon, which has now served its purpose
    loc.likelihood <- loc.likelihood[ -c(1,nrow(loc.likelihood)), c(2,1, 3), drop = FALSE ] ##remove both of the dummy exons, which have now served their purpose

  ################################ Now make it look better, add relevant info
    if (nrow(my.calls$calls) > 0) {

      my.calls$calls$start <- loc.annotations$start[ my.calls$calls$start.p ]
      my.calls$calls$end <- loc.annotations$end[ my.calls$calls$end.p ]
      my.calls$calls$chromosome <- as.character(loc.annotations$chromosome[ my.calls$calls$start.p ])

        my.calls$calls$id <- paste('chr', my.calls$calls$chromosome, ':',  my.calls$calls$start, '-',  my.calls$calls$end, sep = '')
        my.calls$calls$id <- gsub(pattern = "chrchr", replacement = "chr", my.calls$calls$id) ## in case the user already included chr
        my.calls$calls$type <- c('deletion', 'duplication')[ my.calls$calls$type ]

########## make things pretty
      my.calls$calls$BF <- NA
      my.calls$calls$reads.expected <- NA
      my.calls$calls$reads.observed <- NA


      for (ir in 1:nrow(my.calls$calls)) {

        if (my.calls$calls$type[ir] == 'duplication') my.calls$calls$BF[ir] <-  sum(loc.likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir],3 ] - loc.likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir],2 ])

        if (my.calls$calls$type[ir] == 'deletion') my.calls$calls$BF[ir] <-  sum(loc.likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir], 1 ] - loc.likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir],2  ])

        my.calls$calls$reads.expected[ ir ] <-  sum( loc.total [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] * loc.expected [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ ir ] ])
        my.calls$calls$reads.observed[ ir ] <-  sum( loc.test [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] )
      }

      my.calls$calls$reads.expected <- as.integer( my.calls$calls$reads.expected)
      my.calls$calls$reads.ratio <-  signif(my.calls$calls$reads.observed / my.calls$calls$reads.expected, 3)
      my.calls$calls$BF <- signif( log10(exp(1))*my.calls$calls$BF, 3)

      #### shift the numbering properly
      my.calls$calls$start.p <- my.calls$calls$start.p + shift
      my.calls$calls$end.p <- my.calls$calls$end.p + shift

      if (nrow(final) == 0) {final <- my.calls$calls} else {final <- rbind.data.frame(final, my.calls$calls)}
      message('Number of calls for chromosome ', chrom, ' : ', nrow(my.calls$calls))
    }
    shift <- shift + length(good.pos)

  }
  x@CNV.calls <- final
  return (x)
})



#' @title somatic.CNV.call
#' @description Call somatic variants between healthy and disease tissues.
#'
#' @details Use read depth data from targeted sequencing experiments to call CNV between
#' a tumor and matched healthy tissue. This is an experimental function at this stage.
#'
#' @param normal Read count data (numeric vector) for the normal tissue.
#' @param tumor Read count data (numeric vector) for the tumor.
#' @param prop.tumor Proportion of the tumour DNA in the tumour sample (between 0 and 1, and less than 1 if there is normal tissue in the tumor sample).
#' @param chromosome Chromosome information for the bins.
#' @param start Start position of each bin (typically in bp).
#' @param end End position of each bin.
#' @param names Names for each bin (tyically exon names but any way to track
#' the bins will do).
#' @return An ExomeDepth object with CNV calls.
#' @note Absolutely experimental, not the main function from the package.



somatic.CNV.call <- function(normal, tumor, prop.tumor = 1, chromosome, start, end, names) {

  message('Warning: this function is largely untested and experimental')
  message('Initializing the exomeDepth object')
  myTest <- new('ExomeDepth',
                test= tumor,
                reference = normal,
                prop.tumor = prop.tumor,
                formula = 'cbind(test, reference) ~ 1')

  message('Now calling the CNVs')
  myTest <- CallCNVs(x = myTest,
                     transition.probability = 10^-4,
                     chromosome = chromosome,
                     start = start,
                     end = end,
                     name = names)

 return (myTest)
}





#' @title get_loglike_matrix
#' @name  get_loglike_matrix
#' @description Computes the loglikelihood matrix for the three states and each exon
#' @return A likelihood matrix with the states as rows and the exons as columns

NULL



#' @title C_hmm
#' @name C_hmm
#' @description Implements the hidden Markov model (Viterbi algorithm) using a C routine
#' @return A list with two objects: the first contains the optimum Viterbi path and the second the actual CNV calls

NULL
