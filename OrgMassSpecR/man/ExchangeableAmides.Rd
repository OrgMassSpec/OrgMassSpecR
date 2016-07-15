\name{ExchangeableAmides}
\alias{ExchangeableAmides}
\title{Determine the Number of Backbone Amide Hydrogens}
\description{Determine the number of backbone amide hydrogens given a protein/peptide sequence.  Used in hydrogen-deuterium exchange experiments.
}
\usage{ExchangeableAmides(sequence)}
\arguments{
    \item{sequence}{character vector containing one or more amino acid sequences.}
}
\details{
The number of backbone amide hydrogens in an amino acid sequence is the number of residues, minus the number of prolines, minus 1.
}
\value{A numeric vector containing the number of exchangeable hydrogens.}
\author{Nathan G. Dodder}
\seealso{\code{\link{IsotopicDistributionHDX}}}
\examples{
ExchangeableAmides(c("VDVMCTA", "VSTPTL"))

## find the number of exchangeable amides for 
## each peptide in a digest of human serum albumin
x <- Digest(example.sequence, enzyme = "pepsin", IAA = FALSE, missed = 4)
transform(x, exchange = ExchangeableAmides(x$peptide))
}
