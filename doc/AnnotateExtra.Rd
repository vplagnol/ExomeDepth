\name{AnnotateExtra}
\alias{AnnotateExtra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Add annotations to a ExomeDepth object
}
\description{
  Takes annotations in the GRanges format and adds these to the CNV
  calls in the ExomeDepth object.
}
\usage{
AnnotateExtra(x, reference.annotation, min.overlap = 0.5, column.name = "overlap")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An ExomeDepth object.
}
  \item{reference.annotation}{
    The list of reference annotations in GRanges format.
}
  \item{min.overlap}{
    The minimum fraction of the CNV call that is covered by the
    reference call to declare that there is a significant overlap.
}
  \item{column.name}{
    The name of the column used to store the overlap (in the slot CNV.calls).
}
}
\value{
  A ExomeDepth object with the relevant annotations added to the
  CNVcalls slot.
}
\details{A recent version of GenomicRanges (> 1.8.10) is required. Otherwise the
  function will return a warning and not update the ExomeDepth object.}
\references{
  A robust model for read count data in exome sequencing experiments and
  implications for copy number variant calling, Plagnol et al, Bioinformatics
}
\author{
  Vincent Plagnol
}
