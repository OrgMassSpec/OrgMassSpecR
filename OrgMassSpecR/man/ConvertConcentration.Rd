\name{ConvertConcentration}
\alias{ConvertConcentration}
\title{Convert Concentration Basis}
\description{Change the unit basis for a sample concentration, such as
  ng/g wet weight to ng/g dry weight, or pg/g lipid weight to
  pg/g wet weight.
}
\usage{ConvertConcentration(x, convert, percent)
}
\arguments{
  \item{x}{numeric value. The concentration to be converted.}
  \item{convert}{character string. One of \code{wet.to.dry},
    \code{dry.to.wet}, \code{wet.to.lipid}, or \code{lipid.to.wet}. See
    details.}
  \item{percent}{numeric value. Either the percent moisture or percent
    lipid. See details.}
}
\value{Numeric vector of length 1. The unit prefixes are the same as for
  \code{x}, only the concentration basis is converted. 
}
\details{The \code{convert} argument specifies the type of conversion,
  for example \code{wet.to.dry} is wet weight basis to dry weight basis.

  Conversion between wet weight basis and dry weight basis requires
  \code{percent} to equal the percent moisture in the sample. Conversion
  between wet weight basis and lipid weight basis requires
  \code{percent} to equal the percent lipid in the sample. The following
  definitions for percent moisture and percent lipid are used.

  \deqn{percent moisture = \frac{wet weight - dry weight}{wet weight} *
  100}{percent moisture = ((wet weight - dry weight) / wet weight) * 100}

  \deqn{percent lipid = \frac{lipid weight}{wet weight} * 100}{percent
  lipid = (lipid weight / wet weight) * 100} 

  The unit prefixes (for example ng/g or pg/g) are not changed by the
  conversion.
  
}
\examples{
## Convert a concentration of 15.3 ng/g wet weight,
## 5 percent lipid, to ng/g lipid weight.

ConvertConcentration(15.3, "wet.to.lipid", 5)
}
\author{Nathan G. Dodder}

