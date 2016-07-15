\name{FragmentPeptide}
\alias{FragmentPeptide}
\title{Predict Peptide Fragment Ions}
\description{Determine the b- and y-ions or c- and z-ions produced by the fragmentation of a peptide by tandem mass spectrometry.
}
\usage{
FragmentPeptide(sequence, fragments = "by", IAA = TRUE, 
                N15 = FALSE, custom = list())
}
\arguments{
    \item{sequence}{a vector of character strings representing the amino acid sequences to be fragmented by the mass spectrometer.}
    \item{fragments}{character string specifying the fragmentation rules.  Options are \code{"by"} (default) for the b- and y-ions, or \code{"cz"} for the c- and z-ions.}
    \item{IAA}{logical. \code{TRUE} specifies iodoacetylated cysteine and \code{FALSE} specifies unmodified cystine.}
    \item{N15}{logical indicating if the nitrogen-15 isotope should be used in place of the default nitrogen-14 isotope.}
    \item{custom}{a list specifying user defined residues as \code{custom = list(code, mass)}, where \code{code} is a vector of one letter characters and \code{mass} is a vector of the respective monoisotopic masses.  See Details and Examples.}

}
\details{
The amino acid residues must be specified by the one letter codes defined in the help for \code{\link{Digest}}. 

The fragmentation rules can be set for collision induced dissociation (b- and y-ions) or electron transfer dissociation (c- and z-ions). 

The argument \code{IAA} specifies treatment of the protein with iodoacetamide to break the disulfide bonds.  This treatment produces iodoacetylated cysteine residues (elemental formula C5H8N2O2S).

If TRUE, the argument \code{N15} specifies 100\% nitrogen-15 incorporation.  It is intended for proteins grown with a nitrogen-15 labeled food source.  (Although the experiment itself may grow a protein with less than 100\% nitrogen-15 incorporation).  Setting \code{N15 = TRUE} does not modify the mass of a custom residue, or the mass of the nitrogen(s) added if \code{IAA = TRUE}.

If a custom residue code is identical to a predefined residue code, the custom residue mass will be used in place of the predefined mass.
}
\references{
The relative atomic masses of the isotopes are from the NIST Physical Reference Data Website \url{http://physics.nist.gov/PhysRefData/Compositions/}.  The molar mass of a proton (H+) is from the NIST CODATA Website \url{http://physics.nist.gov/cgi-bin/cuu/Value?mmp}.
}
\value{
A data frame with the following column names.  The data frame is arranged this way to facilitate selection of product-precursor ion pairs. 

    \item{ms1seq}{precursor ion sequence.}
    \item{ms1z1}{monoisotopic \emph{m/z} value for the \eqn{[M + H]^{1+}}{[M + H]1+} precursor ions (where M is the precursor mass).}
    \item{ms1z2}{monoisotopic \emph{m/z} value for the \eqn{[M + 2H]^{2+}}{[M + 2H]2+} precursor ions.}
    \item{ms1z3}{monoisotopic \emph{m/z} value for the \eqn{[M + 3H]^{3+}}{[M + 3H]3+} precursor ions.}
    \item{ms2seq}{product ion sequence.}
    \item{ms2type}{the type and charge state of the product ions.}
    \item{ms2mz}{monoisotopic \emph{m/z} values for the product ions.}
}
\author{Nathan G. Dodder and Katharine M. Mullen}
\seealso{\code{\link{MonoisotopicMass}}, \code{\link{Digest}}, \code{\link{PeptideSpectrum}}}
\examples{
## fragment unlabeled peptide
FragmentPeptide("NECFLQHK")

## fragment peptide with carbon-13 labeled lysine
k.mass <- MonoisotopicMass(formula = list(C = 6, H = 12, N = 2, O = 1),
                       isotopes = list(C = 13.0033548378)) 
FragmentPeptide("NECFLQHk", custom = list(code = "k", mass = k.mass))

## fragment peptide with two modifications
m.mass <- MonoisotopicMass(formula = list(C=5, H=9, N=1, O=2, S=1))
FragmentPeptide("NDmELWk", custom = list(code = c("m", "k"), mass = c(m.mass, k.mass)))

## fragment a vector of peptides produced by Digest
x <- Digest(example.sequence)
y <- subset(x, nchar(x$peptide) > 5 & nchar(x$peptide) < 12)
FragmentPeptide(y$peptide)
}
