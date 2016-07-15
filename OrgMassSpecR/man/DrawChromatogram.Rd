\name{DrawChromatogram}
\alias{DrawChromatogram}
\title{Plot a Chromatogram}
\description{Plot a chromatogram, color the area under specified peak(s), and calculate the peak area(s).}
\usage{
DrawChromatogram(time, intensity, range = list(start, stop), color = "blue", 
                 xlab = "retention time", ylab = "intensity", 
                 ylim = c(0, max(intensity) * 1.1), las = 1, ...)
}
\arguments{
    \item{time}{numeric vector containing the time points (the x-axis).}
    \item{intensity}{numeric vector containing the respective signal intensities at each time point (the y-axis).}
    \item{range}{list describing the start and stop time points for each peak, defined as \code{range = list(start, stop)}, where \code{start} is the numeric vector of starting time points and \code{stop} is the numeric vector of the respective ending time points.}
    \item{color}{vector of character strings specifying the color for each peak given in \code{range}.  If a single color is specified, it is applied to all peaks.}
    \item{xlab}{character string specifying the x-axis label.}
    \item{ylab}{character string specifying the y-axis label.}
    \item{ylim}{numeric vector of length 2 specifying the range of the y axis.}
    \item{las}{numeric value specifying the rotation of the axis labels, see par for options.}
    \item{...}{additional parameters to be passed to plot().}
}
\details{The area under the peak(s) is rendered using the \code{\link{polygon}} function.  The area calculation assumes that the polygon does not self-intersect.}
\value{
A data frame with the following column names.
    \item{retention.time}{retention times of the peaks specified in \code{range}}
    \item{area}{the respective peak areas}
    \item{apex.intensity}{the respective intensities at each peak apex} 
}
\author{Nathan G. Dodder}
\examples{
## single peak
x <- DrawChromatogram(example.chromatogram.single$time,  
                      example.chromatogram.single$intensity/100, 
                      range = list(start = 25.4, stop = 26.1), 
                      main = "example chromatogram 1",
                      ylab = "intensity x 100")
# label peak with retention time and area
text(x$retentionTime, x$apexIntensity + 1500, 
     labels = paste("RT = ", round(x$retentionTime, digits = 1),
                    ", Area = ", round(x$peakArea), sep = ""), cex = 0.9) 

## multiple peaks 
y <- DrawChromatogram(example.chromatogram.multiple$time, 
                      example.chromatogram.multiple$intensity / 1000, 
                      range = list(start = c(21.5, 21.925, 23.1, 25.5, 27.35), 
                        stop = c(21.925, 22.4, 23.6, 26.2, 28.0)),
                      color = c("blue", "red", "green", "yellow", "orange"), 
                      xlab = "retention time (min)", 
                      ylab = "intensity x 1000 (cps)", 
                      main = "Example Chromatogram")
## label peaks
text(y$retentionTime, y$apexIntensity + 50, labels = c("a", "b", "c", "d", "e"))
}

