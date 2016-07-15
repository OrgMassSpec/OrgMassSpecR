\name{WriteMspFile}
\alias{WriteMspFile}
\title{Write Spectra to the NIST MSP Text Format}
\description{Writes all spectra in the OrgMassSpecR data frame format to
  a single text file in NIST MSP format. The MSP files can be imported
  into the NIST MS Search program to make a custom seachable library. 
}
\usage{WriteMspFile(spectra, metadata, filename = "library.msp", comment = "")
}
\arguments{
  \item{spectra}{data frame containing the mass spectra.}
  \item{metadata}{data frame containing the metadata.}
  \item{filename}{character string. Name (including path if needed) of
    the MSP text file.}
  \item{comment}{charater string applied to the COMMENT field of the MSP
    file. The comment is applied to all spectra.}
}
\details{The spectra and metadata data frame format is
  specified in the help files for the \code{LibraryReport} functions
  in the spectral library packages, for example see
  \code{SpecLibExample::LibraryReport}. The spectra data frame must
  contain at least columns: \code{filename}, \code{mz}, and
  \code{intensity}. The metadata data frame must contain at least columns: \code{filename} and \code{compound}.

  The NAME field for each spectrum in the MSP file is assigned
  based on \code{compound} in the metadata data frame. 

  Newlines are CR+LF since the NIST MS Search program is
  Windows based.
}
\value{An MSP file is written using a \code{\link{file}} connection.
}
\author{Nathan G. Dodder}
\examples{

## Example using SpecLibExample package
\dontrun{
library(SpecLibExample)
WriteMspFile(spectra = spec,
             metadata = meta,
             comment = "example output")

}
}

