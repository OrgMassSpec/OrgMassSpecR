\name{RetentionIndex}
\alias{RetentionIndex}
\title{Non-Isothermal Retention Index Calculation}
\description{Determine the non-isothermal (temperature-programmed) gas chromatographic retention index of a target compound based on the retention times of adjacent n-alkanes.}
\usage{RetentionIndex(n, target, preceding, following)}
\arguments{
  \item{n}{the number of carbon atoms in the n-alkane eluting prior to the target compound.}
  \item{target}{the retention time of the target compound}
  \item{preceding}{the retention time of the n-alkane eluting prior to the target compound.}
  \item{following}{the retention time of the n-alkane eluting after the target compound.}
}
\details{The retention time units can be either seconds or minutes, but should be consistent. The n-alkanes eluting before and after the target compound should differ by one carbon atom.  This equation does not apply to isothermal separations or when the n-alkane series contains only odd or even carbon atoms.}
\references{Van Den Dool H, Kratz PD, \emph{Journal of Chromatography}, 1963, 11, 463-471.}
\value{A vector giving the retention index.}
\author{Nathan G. Dodder}
\examples{RetentionIndex(8, 30, 28, 33)}
