\name{ConvertPeptide}
\alias{ConvertPeptide}
\title{Convert peptide sequence.}
\description{Convert single amino acid codes to an elemental formula or three letter codes.}
\usage{ConvertPeptide(sequence, output = "elements", IAA = TRUE)}
\arguments{
    \item{sequence}{a character string representing the amino acid sequence.}
    \item{output}{a character string specifying the output.  Options are \code{"elements"} (default) for the elemental formula or \code{"3letter"} for the three letter codes.}
    \item{IAA}{logical. \code{TRUE} specifies iodoacetylated cysteine and \code{FALSE} specifies unmodified cystine.  Used only in determining the elemental formula, not the three letter codes.}
}
\details{
The amino acid residues must be specified by the one letter codes defined in the help for \code{\link{Digest}}.

The argument \code{IAA} specifies treatment of the protein with iodoacetamide to break the disulfide bonds.  This treatment produces iodoacetylated cysteine residues (elemental formula C5H8N2O2S).
}
\value{
If \code{output = "elements"}, the value is a list specifying the number of each element.  This list can be used as input to other functions, see the examples.  If output = \code{"3letter"}, the value is a vector of length 1, containing the amino acid sequence in three letter codes.
}
\author{Nathan G. Dodder}
\examples{
ConvertPeptide("SEQENCE", output = "3letter")

# as input to MonoisotopicMass
MonoisotopicMass(formula = ConvertPeptide("SEQENCE"), charge = 1)

# as input to MolecularWeight
MolecularWeight(formula = ConvertPeptide("SEQENCE"))
}
