\name{ListFormula}
\alias{ListFormula}
\title{Convert an Elemental Formula to a List}
\description{Convert a character string representing an elemental
formula to a list representing the elemental formula.  The list can be
used as input to other functions.}
\usage{
ListFormula(elemental.formula)
}
\arguments{
  \item{elemental.formula}{character string representing the elemental formula.}
}
\details{
To maintain compatibility with \code{\link{MolecularWeight}},
\code{\link{MonoisotopicMass}}, and \code{\link{IsotopicDistribution}}
the elemental formula can contain only C,H,N,O,S,P,Br,Cl,F,Si, and
Sn (Sn is a default element in \code{MolecularWeight} only). Elements not in this set will be ignored and a warning will be
generated. The function can handle repeated elements, but not element
multiplication indicated by parenthesis.

Known issue: Lower case letters after numbers, such as the \dQuote{o} in
\dQuote{C12oBr5H5} will not be caught with a warning and the resulting
list will contain the incorrect number of elements.

}
\seealso{\code{\link{ConvertPeptide}}}
\author{Nathan G. Dodder}
\examples{
ListFormula("C14H8Cl4")
ListFormula("C6H5OH")
}

