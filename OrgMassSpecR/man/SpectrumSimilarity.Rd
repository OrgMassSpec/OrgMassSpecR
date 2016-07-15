\name{SpectrumSimilarity}
\alias{SpectrumSimilarity}
\title{Similarity Between Two Mass Spectra}
\description{Generate a head-to-tail plot of the two mass spectra and calculate a similarity score.}
\usage{
SpectrumSimilarity(spec.top, spec.bottom, t = 0.25, b = 10, 
                   top.label = NULL, bottom.label = NULL, 
                   xlim = c(50, 1200))
}
\arguments{
    \item{spec.top}{data frame containing the experimental spectrum's peak list with the \emph{m/z} values in the first column and corresponding intensities in the second}
    \item{spec.bottom}{data frame containing the reference spectrum's peak list with the \emph{m/z} values in the first column and corresponding intensities in the second}
    \item{t}{numeric value specifying the tolerance used to align the \emph{m/z} values of the two spectra.}
    \item{b}{numeric value specifying the baseline threshold for peak identification.  Expressed as a percent of the maximum intensity.}
    \item{top.label}{character string to label the top spectrum.}
    \item{bottom.label}{character string to label the bottom spectrum.}
    \item{xlim}{numeric vector of length 2, defining the beginning and ending values of the x-axis.}
}
\details{
The mass spectral similarity score is calculated as (where \eqn{\cdot}{\%*\%} is the dot product) \deqn{\cos \theta = \frac{u \cdot v}{\sqrt{u \cdot u} \sqrt{v \cdot v}}}{cos \theta = (u \%*\% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))} where \eqn{u} and \eqn{v} are the aligned intensity vectors of the two spectra.  The \code{t} argument is used to align the intensities.  The bottom spectum is used as the reference spectrum, and the \emph{m/z} values of peaks in the top spectrum that are within \code{t} of a reference \emph{m/z} value are paired with that reference peak.  Ideally, a single peak from the top spectrum should be paired with a single peak from the reference spectrum.  Peaks without a match are paired with an intensity of zero.  A data frame showing the aligned \emph{m/z} values, peak intensities from the top spectrum, and peak intensities from the bottom spectrum, is printed in order to check the results.  The unaligned spectra are shown in the plot.

Note that, although both are based on the cosine of the two intensity vectors, the spectral similarity score given by \code{SpectrumSimilarity} is not the same as that given by the NIST MS Search program, described in the reference below. 
}
\references{
"Optimization and Testing of Mass Spectral Library Search Algorithms for Compound Identification," Stein SE, Scott DR, \emph{Journal of the American Society for Mass Spectrometry}, 1994, 5, 859-866.
}
\value{A vector containing the similarity score. The data frame showing the peak alignment is printed. The head-to-tail plot of the two mass spectra is printed.}
\author{Nathan G. Dodder}
\seealso{\code{\link{PeptideSpectrum}}}
\examples{
SpectrumSimilarity(example.spectrum.unknown, example.spectrum.authentic,
                   top.label = "unknown, electron impact", 
                   bottom.label = "derivatized alanine, electron impact",
                   xlim = c(25, 350))
## label peaks
plot.window(xlim = c(25,350), ylim = c(-125, 125))
text(c(73, 147, 158, 232, 260), c(100, 23, 44, 22, 15) + 10,
     c(73, 147, 158, 232, 260), cex = 0.75)
text(c(73, 147, 158, 232, 260), -c(100, 47, 74, 33, 20) - 10, 
     c(73, 147, 158, 232, 260), cex = 0.75)
mtext("Spectrum similarity", line = 1)
}
