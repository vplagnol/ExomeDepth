\name{qbetabinom.ab}
\alias{qbetabinom.ab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile function for the beta-binomial distribution
}

\description{
  Standard qbetabinomial.ab function which is missing from the VGAM package.
}
\usage{
qbetabinom.ab(p, size, shape1, shape2)
}
\arguments{
  \item{p}{Mean value of the beta-binomial distribution.}
  \item{size}{Size of the beta-binomial.}
  \item{shape1}{First parameter of the beta distribution for p.}
  \item{shape2}{Second parameter of the beta distribution for p.}
}
\value{
  A quantile of the distribution.
}
\seealso{
  \code{VGAM} package.
}
