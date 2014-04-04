\name{CallCNVs}
\alias{CallCNVs}
\title{
  Call CNV data from an ExomeDepth object.
}
\description{
  The function must be called on an ExomeDepth object. Likelihood data
  must have been precomputed which should have been done by default when
  the ExomeDepth object was created.
}
\usage{
CallCNVs(x, chromosome, start, end, name, transition.probability = 0.0001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An \code{ExomeDepth} object}
  \item{transition.probability}{
    Transition probability of the hidden Markov Chain from the normal
    copy number state to either a deletion or a duplication. The default
    (0.0001) expect approximately 20 CNVs genome-wide.
}
  \item{chromosome}{
    Chromosome information for each exon (factor).
}
  \item{start}{
    Start (physical position) of each exon (numeric, must have the same
    length as the chromosome argument).
}
  \item{end}{
    End (physical position) of each exome (numeric, must have the same
    length as the chromosome argument).
}
  \item{name}{
    Name of each exon (character or factor).
  }
}
\details{
  This function fits a hidden Markov model to the read depth data with
  three hidden states (normal, deletion, duplication).
}
\value{
  A list with the following elements:
  \item{calls }{The list of CNV calls, including Bayes factor in favour
    of the call, expected and observed number of reads and position of the
    CNV call.}
}
\author{
  Vincent Plagnol
}
