\name{ReadMspDirectory}
\alias{ReadMspDirectory}
\title{Read in MSP Formatted Mass Spectra}
\description{Reads in all \code{.msp} files within a directory and makes a single concatenated data frame of \emph{m/z} values and intensities.
}
\usage{
ReadMspDirectory(directory, skip = 2, comment.char = "")
}
\arguments{
  \item{directory}{character string. The path to the directory
    containing the \code{.msp} files. Can be \dQuote{./} for the current
    directory on unix-like systems.}
  \item{skip}{integer. The number of lines at the top of the input file to skip before reading. Passed internally to read.table.}
  \item{comment.char}{a character vector of length one containing a single character or an empty string. Use "" to turn off the interpretation of comments altogether. Passed internally to read.table.}
}
\details{MSP is a National Institute of Standards and Technology (NIST) text format for centroid spectra.  This format is exported by some instrument software.

The spectral files do not need to have the extension \code{.msp}. But an attempt will be made to process all files in the directory, regardless of the extension.  The presense of other file types may cause errors.

\code{ReadMspDirectory} is a wrapper for \code{\link{ReadMspFile}}.
}
\value{
A data frame with the following column names:
  \item{mz}{the centroid \emph{m/z} values.}
  \item{intensity}{the respective intensities.}
  \item{filename}{the corresponding filename (used to identify the spectra within the data frame).}
}
\author{Nathan G. Dodder}
\seealso{\code{\link{ReadMspFile}}}
\examples{
## The package directory msp contains two .msp files 
ReadMspDirectory(paste(system.file(package = "OrgMassSpecR"), 
  "/extdata/msp", sep = ""))
}
