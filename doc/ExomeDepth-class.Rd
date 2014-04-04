\name{ExomeDepth-class}
\Rdversion{1.1}
\docType{class}
\alias{ExomeDepth-class}
\alias{AddAnnotations,ExomeDepth-method}
\alias{CallCNVs,ExomeDepth-method}
\alias{TestCNV,ExomeDepth-method}

\title{Class \code{ExomeDepth}}
\description{
  A class to hold the read count data that is used by ExomeDepth to call CNVs.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ExomeDepth",
    data = NULL, test, reference, formula = 'cbind(test, reference) ~
    1', subset.for.speed = NULL)}.
  \code{data} is optional and is only used if the \code{formula}
  argument refers to covariates (in which case these covariates must be
  included in the data frame).
  \code{test} and \code{reference} refer to the read count data for the
  test and reference samples.
  Creating a ExomeDepth object will automatically fit the beta-binomial
  model (using routines from the \code{aod} package) and compute the
  likelihood for the three copy number states (normal, deletion and
  duplication). 
}
\section{Slots}{
  \describe{
    \item{\code{test}:}{\code{numeric}, read count data
      for the test sample.}
    \item{\code{reference}:}{\code{numeric}, read count data for the
      reference sample (usually a combination of samples).} 
    \item{\code{formula}:}{\code{character}, a character string
      describing the linear model linking test and reference. Typically
      this would be \code{cbind(test, reference)  ~ 1}.}
    \item{\code{expected}:}{The expected read count data for the test
      sample assuming normal copy number. }
    \item{\code{phi}:}{The over-dispersion parameter of the binomial
      model. See the \code{aod} package for more details.}
    \item{\code{likelihood}:}{A matrix of likelihood values, one column
      per copy number (deletion, normal, duplication).}
    \item{\code{annotations}:}{A \code{data.frame} specifying the
      chromosome, start and end for the bins used in the read count
      computation.}
    \item{\code{CNV.calls}:}{A \code{data.frame} describing
      the output of the CNV calling procedure.}
  }
}
\section{Methods}{
  \describe{
    \item{CallCNVs}{\code{signature(x = "ExomeDepth",
	transition.probability = "numeric", chromosome = "factor", start
	= "numeric", end = "numeric", name = "character")
      )}: Uses the pre-computed likelihood values and fits a hidden
    Markov Chain to the data to generated merged CNV calls.}
    \item{AddAnnotations}{\code{signature(object = "ExomeDepth", name =
	"character", chromosome = "factor", start = "numeric", end =
	"numeric")}: 
      This method is unlikely to be directly used but it can include the
      exon names, chromosome, start, end into the ExomeDepth object. }
    \item{TestCNV}{\code{signature(x = "ExomeDepth", chromosome =
	"factor", start = "numeric", end = "numeric", type =
	"character")}: type must be either deletion of duplication. This
    function takes an ExomeDepth object and returns the Bayes factor in
    favor of a CNV at the specified location. }
	 }
}
\references{
  Paper recently submitted
}
\author{
  Vincent Plagnol
}
\seealso{
  \code{select.reference.set}
  \code{CallCNVs}
  \code{aod}
}
\examples{
showClass("ExomeDepth")
}
\keyword{classes}
