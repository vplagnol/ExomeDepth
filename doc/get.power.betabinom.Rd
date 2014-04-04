\name{get.power.betabinom}
\alias{get.power.betabinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimate the power to compare two beta-binomial distributions.
}
\description{
  A power study useful in the context of ExomeDepth.
}
\usage{
get.power.betabinom(size,
my.phi,
my.p,
my.alt.p,
theory = FALSE,
frequentist = FALSE,
limit = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{Number of samples from the beta-binomial distribution.}
  \item{my.phi}{Over-dispersion parameter.}
  \item{my.p}{Expected p under the null.}
  \item{my.alt.p}{Expected p under the alternative.}
  \item{theory}{ \code{logical}, should a theoretical limit (large
    sample size) be used? Defaults to FALSE.}
  \item{frequentist}{\code{logical}, should a frequentist version be
    used? Defaults to FALSE.}
  \item{limit}{\code{logical}, should another large sample size limit be
    used? Defaults to FALSE.}
}
\value{
  An expected Bayes factor.
}
\author{
  Vincent Plagnol
}
