\name{IsotopicDistributionHDX}
\alias{IsotopicDistributionHDX}
\title{Isotopic Distribution of a Peptide Undergoing H-D Exchange}
\description{
Simulates the isotopic distribution of a peptide undergoing hydrogen-deuterium exchange.  Only the peptide backbone amides are labeled with deuterium.
}
\usage{
IsotopicDistributionHDX(sequence, incorp, charge = 1, 
                        custom = list(code = NULL, elements = NULL))
}
\arguments{
    \item{sequence}{character vector specifying the amino acid sequence.}
    \item{incorp}{numeric value from 0 to 1, specifying the fraction of deuterium incorporation in the backbone amides.}
    \item{charge}{numeric value specifying the number of positive charges (\eqn{[M + nH]^{n+}}{[M + nH]n+})}
    \item{custom}{a list specifying user defined residues as \code{custom = list(code, mass)}, where \code{code} is a vector of one letter characters and \code{mass} is a vector of the respective monoisotopic masses.}
}
\details{
The amino acid residues must be specified by the one letter codes defined in the help for \code{\link{Digest}}.  If a custom residue code is identical to a predefined residue code, the custom residue mass will be used in place of the predefined mass.

The natural incorporation of deuterium is 0.000115.

The algorithm used in this function is based on \code{\link{sample}}, and will give a slightly different result each time it is run.  \code{IsotopicDistributionHDX} has not been tested for sequences over approximately 3000 amu and charges greater that 3+.
}
\references{The relative atomic masses of the isotopes are from the NIST Physical Reference Data Website \url{http://physics.nist.gov/PhysRefData/Compositions/}.  The molar mass of a proton (H+) is from the NIST CODATA Website \url{http://physics.nist.gov/cuu/Constants/index.html}.  
}
\value{
A data frame with the following column names.
    \item{mz}{the \emph{m/z} value}
    \item{intensity}{the number of counts at each \emph{m/z} value.  The total is 10000.}
    \item{percent}{the intensity at each \emph{m/z} value, expressed as a percent of the maximum intensity.}
}
\author{Nathan G. Dodder}
\seealso{\code{\link{IsotopicDistribution}}, \code{\link{IsotopicDistributionN}}, \code{\link{sample}}}
\examples{
## simulate a peptide with 0.0115 and 100 percent incorporation and plot
x <- IsotopicDistributionHDX("NECFLQHK", incorp = 0.000115)
x$t <- "incorp = 0.0115\%"
y <- IsotopicDistributionHDX("NECFLQHK", incorp = 1)
y$t <- "incorp = 100\%"
z <- rbind(x, y)
library(lattice)
print(xyplot(percent ~ mz | t, data = z,
       type = "h", main = "H-D exchange simulation",
       xlab = "m/z", ylab = "intensity (\%)"))
}
