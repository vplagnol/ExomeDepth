
#' Plotting function for ExomeDepth objects
#'
#' Plot function for the ExomeDepth class
#'
#' @name plot-methods
#' @aliases plot-methods plot.ExomeDepth plot,ANY-method
#' plot,ExomeDepth,ANY-method plot,ExomeDepth-method
#' @docType methods
#' @param x ExomeDepth object
#' @param sequence character, Name of the sequence/chromosome of the region to plot (for example "chr5" would be typical)
#' @param xlim numeric of size 2, start and end position of the region to plot
#' @param ylim numeric of size 2, range for the y-axis
#' @param count.threshold numeric, minimum number of reads in the reference set to display a point in the plot
#' @param ylab Defaults to ''
#' @param xlab Defaults to ''
#' @param type Defaults to 'b'
#' @param pch Defaults to '+'
#' @param with.gene Logical, defaults to FALSE, Should the gene information (obtained from the annotation data) be plotted under the read depth plot?
#' @param col character, Colour for the line displaying the read depth ratio for each exon
#' @param ... Additional arguments to be passed to the base plot function



setMethod("plot",
          "ExomeDepth",
          function(x, sequence, xlim, ylim = NULL, count.threshold = 10, ylab = 'Observed by expected read ratio', xlab = '', type = 'b', pch = '+', with.gene = FALSE, col = 'red', ...) {

  anno <- x@annotations
  selected <-  which(anno$chromosome == sequence & anno$start >= xlim[1] & anno$end <= xlim[2] & (x@test + x@reference)*x@expected > count.threshold)

  if (length(selected) == 0) {
    warning('No exon seem to be located in the requested region. It could also be that the read count is too low for these exons? In any case no graph will be plotted.')
    return(NULL)
  }


### From now on we can assume the selection is not empty
  anno <- anno[selected,]

  anno$expected <- x@expected[ selected ]
  anno$freq <- x@test[ selected ]/ (x@reference[selected ] + x@test[selected])
  anno$middle <- 0.5*(anno$start + anno$end)
  anno$ratio <-  anno$freq/ anno$expected
  anno$test <- x@test[ selected ]
  anno$reference <- x@reference[ selected ]
  anno$total.counts <- anno$test + anno$reference
  if (length( x@phi ) == 1) anno$phi <- x@phi else anno$phi <-  x@phi [ selected ]

  #browser()
  for (i in 1:nrow(anno)) {
    anno$my.min.norm[ i ] <- qbetabinom (p = 0.025, size =  anno$total.counts[ i ], phi = anno$phi[ i ], prob =  anno$expected[ i ])
    anno$my.max.norm[ i ] <- qbetabinom (p = 0.975, size =  anno$total.counts[ i ], phi = anno$phi[ i ], prob =  anno$expected[ i ])
  }

  #anno$my.min.norm.prop <-   anno$my.min.norm / (anno$my.min.norm + anno$reference)
  #anno$my.max.norm.prop <-   anno$my.max.norm / (anno$my.max.norm + anno$reference)

  anno$my.min.norm.prop <-   anno$my.min.norm /  anno$total.counts
  anno$my.max.norm.prop <-   anno$my.max.norm /  anno$total.counts


######### Now we plot
  if (with.gene) {  ##if we want the gene information we add this extra bit
    layout(mat = matrix(data = 1:2, nrow = 2, ncol = 1), widths = c(1, 1), heights = c(2, 1) )
    par(mar = c(0, 4, 2, 2))
  }

  if (is.null(xlim)) xlim <- range(anno$middle)
  if (is.null(ylim)) ylim <- c( 0, max( anno$ratio, 1.25*max( anno$my.max.norm.prop/anno$expected, na.rm = TRUE )))

  plot (x = NA,
        y = NA,
        ylim = ylim,
        xlim = xlim,
        ylab = ylab,
        xaxt = ifelse (with.gene, 'n', 's'),
        ...)

  #title(main = 'A', adj = 0)
  polygon(x = c(anno$middle , rev(anno$middle)),
          y = c(anno$my.min.norm.prop, rev(anno$my.max.norm.prop))/c(anno$expected, rev(anno$expected)),
          col = 'grey',
          border = NA)

  ###### make the borders pretty
  polygon(x = c(xlim[1] , anno$middle[1], anno$middle[1], xlim[1]),
          y = c(rep(anno$my.min.norm.prop[1], 2), rep(anno$my.max.norm.prop[1], 2))/anno$expected[1],
          col = 'grey',
          border = 'grey')

  nr <- nrow(anno)
  polygon(x = c(xlim[2] , anno$middle[nr], anno$middle[nr], xlim[2]),
          y = c(rep(anno$my.min.norm.prop[nr], 2), rep(anno$my.max.norm.prop[nr], 2))/anno$expected[nr],
          col = 'grey',
          border = 'grey',
          lty = 1)




  points (x = anno$middle[ order(anno$middle) ],
          y = anno$ratio [ order(anno$middle) ],
          type = 'b',
          pch = '+',
          col = col,
          ...)


  if (with.gene) {

    message("Plotting the gene data")
    par(mar = c(4,4,0,2))

    plot (x = NA,
          y = NA,
          ylim = c(0,1),
          xlim = xlim,
          yaxt = 'n',
          ylab = '',
          xlab = xlab,
          xaxt = 'n')

    ##### This bit to make sure that the axis is plotted with integer values
    my.pos <- axTicks(side = 1)
    axis(side = 1, at = my.pos, labels = as.integer(my.pos))


    ################ Now subset the annotations we want to show


    #if (is.null(annotations)) {
    #anno <- x@annotations
    selected <-  which(anno$chromosome == sequence & anno$start >= xlim[1] & anno$end <= xlim[2])
    selected.2 <- max(min(selected) - 1, 1) :  min( max(selected) + 1, nrow(anno))   ##complicated line that adds one exon on each side of the window
    exon.array <- anno[ selected.2,,drop = FALSE]
   # } else {
   #   selected <-  which(annotations$chromosome == sequence & anno$start >= xlim[1] & annotations$end <= xlim[2])
   #   exon.array <- annotations[ selected.2,,drop = FALSE]
   # }


    exon.array$short.name <- gsub(exon.array$name, pattern = '-.*', replacement = '')

### Now if an additional gene was added at either end, remove it
      if (nrow(exon.array) > 1) {
          if (exon.array$short.name[1] != exon.array$short.name[2]) {
              exon.array <- exon.array[-1, ]
          }
          if (nrow(exon.array)>1 && (exon.array$short.name[nrow(exon.array)] !=
                                     exon.array$short.name[nrow(exon.array) - 1])) {
              exon.array <- exon.array[-nrow(exon.array), ]
          }
      }

    exon.array$start.gene <- tapply(INDEX = exon.array$short.name, exon.array$start, FUN = min) [ exon.array$short.name ]
    exon.array$middle <- 0.5*(exon.array$start + exon.array$end)
    exon.array <- exon.array[ exon.array$short.name != 'RP11' &
                             ! grepl(pattern = 'ENST.*', exon.array$short.name) &
                             ! grepl(pattern = 'AC0.*', exon.array$short.name) ,]



    if (nrow(exon.array) >= 1) {
      pos <- 1
      lev <- 0.4
      arr <- by (exon.array, exon.array$start.gene, FUN = function(x){
        my.x <- range(x$middle)
        lines(x = my.x, y = c(lev, lev))

        my.x <- mean(range(x$middle))
        if (my.x < xlim[1]) my.x <- xlim[1]
        if (my.x > xlim[2]) my.x <- xlim[2]

        text(x = my.x,
             y = ifelse (lev == 0.6, 0.7, 0.3),
             labels = x$short.name[1],
             pos = pos,
             cex = 0.5)


        for (i in 1:nrow(x)) {
          lines(x = c(x$middle[i], x$middle[i]), y = c(lev - 0.1, lev + 0.1))
        }

        if (pos == 1) {pos <<- 3; lev <<- 0.6;} else {pos <<- 1; lev <<- 0.4}  ##this is for the text

      })
    }
  }
})
