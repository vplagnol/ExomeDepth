

setGeneric("AnnotateExtra", def = function(x, reference.annotation, min.overlap = 0.5, column.name) standardGeneric('AnnotateExtra'))


#' @title AnnotateExtra
#' @aliases AnnotateExtra
#' @description Add annotations to a ExomeDepth object.
#' @md
#' @details This function takes annotations in the GRanges format and adds these to the CNV calls in
#' the ExomeDepth object.
#' Note that a recent version of `GenomicRanges (> 1.8.10)` is required. Otherwise the
#' function will return a warning and not update the ExomeDepth object.
#'
#' @param x An `ExomeDepth` object.
#' @param reference.annotation The list of reference annotations in GRanges format.
#' @param min.overlap Numeric, defaults to 0.5. This defines the minimum fraction of the CNV call that is covered by the reference call to declare that there is a significant overlap.
#' @param column.name The name of the column used to store the overlap (in the slot CNV.calls).
#' @return An ExomeDepth object with the relevant annotations added to the CNVcalls slot.

setMethod("AnnotateExtra", "ExomeDepth", function( x, reference.annotation, min.overlap, column.name) {

    my.calls.GRanges <- GenomicRanges::GRanges(seqnames = factor(x@CNV.calls$chromosome),
                                               IRanges::IRanges(start=x@CNV.calls$start,end= x@CNV.calls$end))

    test <- GenomicRanges::findOverlaps(query = my.calls.GRanges,
                                        subject = reference.annotation)
    test <- data.frame(calls = test@from, ref = test@to)

    ##add info about the CNV calls
    test$call.start <- x@CNV.calls$start[ test$calls ]
    test$call.end <- x@CNV.calls$end[ test$calls ]
    test$chromosome.end <- x@CNV.calls$chromosome[ test$calls ]

    ## info about the reference calls
    test$callsref.start <- GenomicRanges::start(reference.annotation) [ test$ref ]
    test$callsref.end<- GenomicRanges::end(reference.annotation) [ test$ref ]

    ## line below purely for CRAN compliance
    callsref.start <- callsref.end <- overlap <- call.end <- call.start <- NULL

    ## estimate the overlap
    test <- dplyr::mutate(test,
                          overlap = pmin (callsref.end, call.end) -  pmax( call.start, callsref.start)) %>%
        dplyr::filter(overlap > min.overlap*(call.end - call.start))

    my.split <-  split(as.character(GenomicRanges::elementMetadata(reference.annotation)$names)[ test$ref], f = test$calls)
    my.overlap.frame <- data.frame(call = names(my.split),  target = sapply(my.split, FUN = paste, collapse = ','))
    my.overlap.frame <- data.frame(call = names(my.split),  target = sapply(my.split, FUN = paste, collapse = ','))


    x@CNV.calls[, column.name] <- as.character(my.overlap.frame$target)[ match(1:nrow(x@CNV.calls), table = my.overlap.frame$call) ]
    return(x)
})

