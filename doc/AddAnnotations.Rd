\name{AddAnnotations}
\alias{AddAnnotations}
\title{
Add annotation data to an exome dataset (i.e. position and names of the exons)
}
\description{
Takes the positions of the exons and the names and include this
information into a ExomeDepth object.
}
\usage{
AddAnnotations(x, name, chromosome, start, end)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An ExomeDepth object.
}
  \item{name}{
    Names of the exons(\code{character}).
}
  \item{chromosome}{
    Chromosome information for each exon (\code{factor}).
}
  \item{start}{
    (numeric) Start of the exons.
}
  \item{end}{
    (numeric) End of the exons.
}
}
\details{
  Unlikely to be used directly but rather as a subroutine of the
  CallCNVs function.
}
\value{
  The same ExomeDepth object with the annotation data frame added.
}
\author{
  Vincent Plagnol
}
\note{
  May not be ran often without a main function call from ExomeDepth. It
  may be useful as s standalone function in some cases however.
}
