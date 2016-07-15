\name{FlowTime}
\alias{FlowTime}
\title{Solvent transit time.}
\description{Calculate the time required for a liquid to flow through a defined length of tubing.  Intended for conventional-flow or nano-flow liquid chromatography.}
\usage{FlowTime(internalDiameterMicrometers = 24, tubingLengthCentimeters = 45, 
                flowRateMicrolitersPerMinute = 0.3)}
\arguments{
    \item{internalDiameterMicrometers}{the internal diameter of the tubing, in micrometers; numeric value.}
    \item{tubingLengthCentimeters}{the length of tubing, in centimeters; numeric value.}
    \item{flowRateMicrolitersPerMinute}{the flow rate of the solvent in microliters per minute; numeric value}
}
\details{The calculation assumes the viscosity of the liquid is negligible.}
\value{Vector of the flow time in seconds.}
\author{Nathan G. Dodder}
\seealso{\code{\link{DeadVolume}}}
\examples{FlowTime(internalDiameterMicrometers = 47, tubingLengthCentimeters = 33, 
                   flowRateMicrolitersPerMinute = 5)}
