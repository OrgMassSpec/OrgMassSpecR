\name{Digest}
\alias{Digest}
\title{Predict Peptides Resulting from Enzymatic Digest}
\description{Cleave an amino acid sequence (a protein or peptide) according to enzyme specific rules and calculate the precursor ion \emph{m/z} values.
}
\usage{
Digest(sequence, enzyme = "trypsin", missed = 0, IAA = TRUE, 
       N15 = FALSE, custom = list())
}
\arguments{
    \item{sequence}{a character string representing the amino acid sequence to be cleaved by the enzyme.}
    \item{enzyme}{a character string specifying the rules for cleavage.  Options are \code{"trypsin"} (default),                    \code{"trypsin.strict"} (see details), or \code{"pepsin"}.}
    \item{missed}{the maximum number of missed cleavages.  Must be an integer of 0 (default) or greater.  An error will result if the specified number of missed cleavages is greater than the maximum possible number of missed cleavages.}
    \item{IAA}{logical. \code{TRUE} specifies iodoacetylated cysteine and \code{FALSE} specifies unmodified cystine.}
    \item{N15}{logical indicating if the nitrogen-15 isotope should be used in place of the default nitrogen-14 isotope.}
    \item{custom}{a list specifying user defined residues as \code{custom = list(code, mass)}, where \code{code} is a vector of one letter characters and \code{mass} is a vector of the respective monoisotopic masses.  See Details and Examples.}
}
\details{The amino acid residues must be specified by one letter codes.  The predefined residues are: 

\tabular{ll}{
A = alanine         \tab L = leucine\cr 
R = arginine        \tab K = lysine\cr
N = asparagine      \tab M = methionine\cr
D = aspartic acid   \tab F = phenylalanine\cr
C = cysteine        \tab P = proline\cr
E = glutamic acid   \tab S = serine\cr
Q = glutamine       \tab T = threonine\cr
G = glycine         \tab W = tryptophan\cr
H = histidine       \tab Y = tyrosine\cr
I = isoleucine      \tab V = valine
}

If \code{"trypsin"} is specified, the sequence is cleaved on the c-terminal side of K and R residues, except if K or R is followed by P.  If \code{"trypsin.strict"} is specified, the sequence is cleaved on the c-terminal side of K and R residues. If \code{"pepsin"} is specified, the sequence is cleaved on the c-terminal side of F, L, W, Y, A, E, and Q residues.  This rule is specific to pepsin at pH \eqn{>} 2, as used in hydrogen-deuterium exchange experiments.

When \code{"trypsin"} is specified, KP and RP are not considered missed cleavages when \code{missed} is \eqn{>} 0.

The argument \code{IAA} specifies treatment of the protein with iodoacetamide.  This treatment produces iodoacetylated cysteine residues (elemental formula C5H8N2O2S).

If \code{TRUE}, the argument \code{N15} specifies 100\% nitrogen-15 incorporation.  It is intended for proteins grown with a nitrogen-15 labeled food source.  (Although the experiment itself may grow a protein with less than 100\% nitrogen-15 incorporation).  Setting \code{N15 = TRUE} does not modify the mass of a custom residue, or the mass of the nitrogen(s) added if \code{IAA = TRUE}.

If a custom residue code is identical to a predefined residue code, the custom residue mass will be used in place of the predefined mass.  

The error message \dQuote{\code{object "mass" not found}} indicates the input sequence contains an undefined residue(s).
}
\references{The relative atomic masses of the isotopes are from the NIST Physical Reference Data Website \url{http://physics.nist.gov/PhysRefData/Compositions/}.  The molar mass of a proton (H+) is from the NIST CODATA Website \url{http://physics.nist.gov/cuu/Constants/index.html}.  
}
\value{A data frame with the following column names.
    \item{peptide}{resulting peptides.}
    \item{start}{beginning residue positions in the the original sequence.}
    \item{end}{ending residue positions in the the original sequence.}
    \item{mc}{number of missed cleavages.}
    \item{mz1}{monoisotopic \emph{m/z} values for the \eqn{[M + H]^{1+}}{[M + H]1+} ions (where M is the precursor mass).}
    \item{mz2}{monoisotopic \emph{m/z} values for the \eqn{[M + 2H]^{2+}}{[M + 2H]2+} ions.}
    \item{mz3}{monoisotopic \emph{m/z} values for the \eqn{[M + 3H]^{3+}}{[M + 3H]3+} ions.}
}  
\author{Nathan G. Dodder}
\seealso{\code{\link{MonoisotopicMass}}, \code{\link{FragmentPeptide}}}
\examples{
## digest human serum albumin with 0 and 1 missed cleavages
Digest(example.sequence, missed = 1)

## digest human serum albumin with a phosphoserine at position 58
## and all methionines oxidized
modifiedHsaSequence <- strsplit(example.sequence, split = "")[[1]]
modifiedHsaSequence[58] <- "s"   # insert code for phosphoserine  
modifiedHsaSequence <- paste(modifiedHsaSequence, collapse = "")
Digest(modifiedHsaSequence, custom = list(code = c("s","M"), 
       mass = c(MonoisotopicMass(list(C=3, H=6, N=1, O=5, P=1)),
                MonoisotopicMass(list(C=5, H=9, N=1, O=2, S=1)))))
                
## digest human serum albumin with strict rules
Digest(example.sequence, enzyme = "trypsin.strict")
}
