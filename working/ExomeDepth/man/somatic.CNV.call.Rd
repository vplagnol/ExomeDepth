\name{somatic.CNV.call}
\alias{somatic.CNV.call}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Call somatic variants between healthy and disease tissues.
}
\description{
  Use read depth data from targeted sequencing experiments to call CNV
  between a tumor and matched healthy tissue.
}
\usage{
somatic.CNV.call(normal, tumor, prop.tumor = 1, chromosome, start, end, names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{normal}{
    Read count data (numeric vector) for the normal tissue.
}
  \item{tumor}{
    Read count data (numeric vector) for the tumor.
}
  \item{prop.tumor}{
    Proportion of the tumour DNA in the tumour sample (between 0 and 1,
    and less than 1 if there is normal tissue in the tumor sample). 
}
  \item{chromosome}{
    Chromosome information for the bins.
}
  \item{start}{
    Start position of each bin (typically in bp).
}
  \item{end}{
    End position of each bin.
}
  \item{names}{
    Names for each bin (tyically exon names but any way to track the
    bins will do).
  }
}
\details{
  Experimental function at this stage.
}
\value{
  An ExomeDepth object with CNV calls.
}
\references{
  A robust model for read count data in exome sequencing experiments and
  implications for copy number variant calling, Bioinformatics, In Press
}
\author{
  V. Plagnol
}
\note{
  Absolutely experimental, not the main function from the package.
}

