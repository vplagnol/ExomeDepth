\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot.ExomeDepth}
\alias{plot,ANY-method}
\alias{plot,ExomeDepth,ANY-method}
\alias{plot,ExomeDepth-method}
\title{Plotting function for ExomeDepth objects}
\description{
Plot function for the ExomeDepth class
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ExomeDepth", sequence, xlim, ylim = NULL,
    count.threshold = 10, ylab = 'Observed by expected read ratio', xlab
    = '', type = 'b', pch = '+', with.gene = FALSE, col =
    'red',...)}}{}}}
\arguments{
  \item{x}{ExomeDepth object}
  \item{sequence}{character, Name of the sequence/chromosome of the region to
    plot (for example "chr5" would be typical)}
  \item{xlim}{numeric of size 2, start and end position of the region
    to plot}
  \item{ylim}{numeric of size 2, range for the y-axis}
  \item{count.threshold}{numeric, minimum number of reads in the reference set to
    display a point in the plot}
  \item{with.gene}{logical, Should the gene information (obtained from
    the annotation data) be plotted under the read depth plot?}
  \item{col}{character, Colour for the line displaying the read depth
    ratio for each exon}
}
\keyword{methods}
