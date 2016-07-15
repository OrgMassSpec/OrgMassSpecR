\name{ReadMspFile}
\alias{ReadMspFile}
\title{Read in a MSP Formatted Mass Spectrum}
\description{Reads in a \code{.msp} file and makes a data frame of \emph{m/z} values and intensities.
}
\usage{
ReadMspFile(file, skip = 2, comment.char = "")
}
\arguments{
  \item{file}{character string. The \code{.msp} filename.  Include the path to the filename if necessary.}
  \item{skip}{integer. The number of lines at the top of the input file to skip before reading. Passed internally to read.table.}
  \item{comment.char}{a character vector of length one containing a single character or an empty string. Use "" to turn off the interpretation of comments altogether. Passed internally to read.table.}
}
\details{MSP is a National Institute of Standards and Technology (NIST) ASCII format for centroid spectra.  This format is exported by some instrument software.

The spectral file does not need to have the extension \code{.msp}.
}
\value{
A data frame with the following column names:
  \item{mz}{the centroid \emph{m/z} values.}
  \item{intensity}{the respective intensities.}
}
\author{Nathan G. Dodder}
\seealso{\code{\link{ReadMspDirectory}}}
\examples{
ReadMspFile(paste(system.file(package = "OrgMassSpecR"), "/extdata/msp/pentaBDE.msp", sep = ""))
}
