\name{TestCNV}
\alias{TestCNV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Computes the Bayes Factor in favour of a CNV defined by position and type.
}
\description{
  Test what evidence supports the presence of a CNV in an ExomeDepth object.
}
\usage{
TestCNV(x, chromosome, start, end, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ExomeDepth object containing the likelihood information.} 
  \item{chromosome}{Chromosome data (\code{factor})}
  \item{start}{Start of the putative CNV}
  \item{end}{End of the putative CNV}
  \item{type}{\code{character}, Should be either 'deletion' or 'duplication'}
}
\value{
  A Bayes factor assessing the evidence in favoor of the CNV.
}
