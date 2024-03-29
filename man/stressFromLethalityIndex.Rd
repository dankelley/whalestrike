% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lethality_index.R
\name{stressFromLethalityIndex}
\alias{stressFromLethalityIndex}
\title{Compute stress, based on lethality index}
\usage{
stressFromLethalityIndex(injury)
}
\arguments{
\item{injury}{numerical value or vector, giving threat of injury (in range 0 to 1).}
}
\value{
whale compression stress, in Pascals.
}
\description{
The model used for this is the logistic model, fitting observed injury/lethality
statistics to the base-10 logarithm of the maximum compression stress during
a simulated impact event.
}
\examples{
stressFromLethalityIndex(0.5) # approx. 254000 Pa, i.e. parameters()$logistic$tau50

}
\seealso{
Other functions dealing with Whale Lethality index: 
\code{\link{lethalityIndexFromStress}()}
}
\author{
Dan Kelley
}
\concept{functions dealing with Whale Lethality index}
