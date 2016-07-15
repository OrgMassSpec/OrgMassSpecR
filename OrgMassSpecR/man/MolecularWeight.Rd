\name{MolecularWeight}
\alias{MolecularWeight}
\title{Calculate the molecular weight of an organic molecule.}
\description{Given an elemental formula and the average relative atomic masses of the elements, determine the molecular weight.}
\usage{MolecularWeight(formula = list(), amu = list())}
\arguments{  
    \item{formula}{a list describing the elemental formula.  The allowed
      elements are C, H, N, O, S, P, Br, Cl, F, Si, Sn, and a user defined "x".  See Examples.}
    \item{amu}{a list specifing user defined standard atomic masses of the elements (in atomic mass units).}
}
\details{
The user defined \code{x} in the argument lists can be used to define an additional element, such as a metal.  See example.

This function will accept values that do not correspond to known physical reality, such as a fractional number of elements or the wrong standard atomic mass for an element.
}
\references{
The relative atomic masses of the elements are from the NIST Physical Reference Data Website \url{http://physics.nist.gov/PhysRefData/Compositions/}.  
}
\value{The molecular weight of the molecule.}
\author{Nathan G. Dodder}
\seealso{\code{\link{MonoisotopicMass}}}
\examples{
MolecularWeight(formula = list(C=2, H=4))

## molecular weight of cyanocobalamin (C63H88CoN14O14P)
MolecularWeight(formula = list(C=63, H=88, N=14, O=14, P=1, x=1),
                amu = list(x = 58.933200))
}
