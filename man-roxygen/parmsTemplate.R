#' @param parms A named list holding model parameters, created by
#' \code{\link{parameters}}.

## \describe{
## \item{ms}{Ship mass [kg].}
## \item{as}{Ship area [m^2]. This, together with \code{CDs}, is used by
#for water drag calculation, using CD=1e-3.}
## \item{B}{Ship impact horizontal extent [m].}
## \item{D}{Ship impact vertical extent [m].}
## \item{mw}{Whale mass [kg]. (Consider using \code{\link{massFromLength}},
##  if length data are easier to obtain than mass data.)}
## \item{delta}{Whale skin thickness [m]. See e.g. Miller et al. (2011), for
##  right whales.}
## \item{Es}{Whale skin elastic modulus [Pa].}
## \item{theta}{Whale skin deformation angle [deg].}
## \item{beta}{Whale blubber thickness [m].}
## \item{blubbermodel}{Blubber model to use (passed to \code{\link{blubberForce}}).}
## }
