\name{qbetabinom}
\alias{qbetabinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Quantile for betabin function
}
\description{
  Quantile function for the betabinomal distribution using the p/phi parameterisation.
}
\usage{
qbetabinom(p, size, phi, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Point of the distribution from which one is looking for the quantile}
  \item{size}{Sample size of the random variable}
  \item{phi}{Over-dispersion parameter}
  \item{prob}{Mean probability of the binomial distribution}
}
\details{
  Filling a gap in the VGAM package.
}
\value{A real number corresponding to the quantile p.}
\author{Vincent Plagnol}
\seealso{
  VGAM R package.
}
