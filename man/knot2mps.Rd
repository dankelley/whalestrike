% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{knot2mps}
\alias{knot2mps}
\title{Convert a speed in knots to a speed in m/s}
\usage{
knot2mps(knot)
}
\arguments{
\item{knot}{Speed in knots.}
}
\value{
Speed in m/s.
}
\description{
See also \code{\link[=mps2knot]{mps2knot()}}, which is the inverse of this function.
}
\examples{
library(whalestrike)
knots <- seq(0, 20)
plot(knots, knot2mps(knots), xlab = "Speed [knots]", ylab = "Speed [m/s]", type = "l")

}
\seealso{
Other functions dealing with units: 
\code{\link{mps2knot}()}
}
\author{
Dan Kelley
}
\concept{functions dealing with units}
