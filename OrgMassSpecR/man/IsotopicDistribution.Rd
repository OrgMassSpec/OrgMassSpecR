\name{IsotopicDistribution}
\alias{IsotopicDistribution}
\title{Isotopic Distribution of an Organic Molecule.}
\description{Simulate the isotopic distribution of an organic molecule using the natural abundances of the isotopes.
}
\usage{
IsotopicDistribution(formula = list(), charge = 1)
}
\arguments{
  \item{formula}{a list describing the charged elemental formula.  The
    allowed elements are C, H, N, O, S, P, Br, Cl, F, Si.}
  \item{charge}{an integer specifying the number of positive or negative
    charges. A charge of zero is not allowed.}
}
\details{
The elemental formula should be that of the charged molecule; i.e, the
\code{charge} argument does not change the elemental formula.

The algorithm used in this function is based on \code{\link{sample}},
and will give a slightly different result each time it is run.

This function is intended for simulating the isotopic distributions of
small molecules (< approximately 3000 amu and <= 3 charges) as measured on low resolution mass spectrometers. 
}
\references{
The relative atomic masses of the isotopes are from the NIST Physical Reference Data Website \url{http://physics.nist.gov/PhysRefData/Compositions/}.
}
\value{
A data frame with the following column names.

    \item{mz}{the \emph{m/z} value}
    \item{intensity}{the number of counts at each \emph{m/z} value.  The total is 10000.}
    \item{percent}{the intensity at each \emph{m/z} value, expressed as a percent of the maximum intensity.}
}
\author{Nathan G. Dodder}
\seealso{\code{\link{IsotopicDistributionN}}, \code{\link{IsotopicDistributionHDX}}, \code{\link{sample}}}
\examples{
x <- IsotopicDistribution(formula = list(C = 9, H = 4, Br = 5, Cl = 1, N = 2))
library(lattice)
xyplot(percent ~ mz, data = x,
  type = "h",
  xlab = "m/z",
  ylab = "intensity (\%)",
  main = "Isotopic Distribution, C9H4Br5ClN2")
}
