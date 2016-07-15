\name{DeadVolume}
\alias{DeadVolume}
\title{Internal volume of tubing.}
\description{Calculate the internal volume of a defined length of tubing.}
\usage{DeadVolume(internalDiameterMicrometers = 24, tubeLengthCentimeters = 45)}
\arguments{
    \item{internalDiameterMicrometers}{the internal diameter of the tubing, in micrometers; numeric value.}
    \item{tubeLengthCentimeters}{the length of tubing, in centimeters; numeric value.}
}
\value{Vector of the dead volume in microliters.}
\author{Nathan G. Dodder}
\seealso{\code{\link{FlowTime}}}
\examples{DeadVolume(internalDiameterMicrometers = 47, tubeLengthCentimeters = 33)}
