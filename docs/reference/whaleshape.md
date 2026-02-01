# Whale side-view shape

This is a data frame containing 45 points specifying the shape of a
right whale, viewed from the side. It was created by digitizing the
whale shape (ignoring fins) that is provided in the necropsy reports of
Daoust et al. (2018). The data frame contains `x` and `y`, which are
distances nondimensionalized by the range in `x`; that is, `x` ranges
from 0 to 1. The point at the front of the whale is designated as x=y=0.

## Usage

``` r
whaleShape()
```

## References

Daoust, Pierre-Yves, Emilie L. Couture, Tonya Wimmer, and Laura Bourque.
"Incident Report. North Atlantic Right Whale Mortality Event in the Gulf
of St. Lawrence, 2017." Canadian Wildlife Health Cooperative, Marine
Animal Response Society, and Fisheries and Oceans Canada, 2018.
<https://publications.gc.ca/site/eng/9.850838/publication.html>.

## Author

Dan Kelley

## Examples

``` r
library(whalestrike)
shape <- whaleShape()
plot(shape$x, shape$y, asp = 1, type = "l")
polygon(shape$x, shape$y, col = "lightgray")
lw <- 13.7
Rmax <- 0.5 * lw * diff(range(shape$y))
mtext(sprintf("Max. radius %.2fm for %.1fm-long whale", Rmax, lw), side = 3)

```
