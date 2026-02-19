# Convert a speed in m/s to a speed in knots

This is done by dividing by the factor 1.852e3/3600, See also
[`knot2mps()`](https://dankelley.github.io/whalestrike/reference/knot2mps.md),
which is the inverse of this function.

## Usage

``` r
mps2knot(mps)
```

## Arguments

- mps:

  Speed in metres per second.

## Value

Speed in knots.

## See also

Other functions dealing with units:
[`knot2mps()`](https://dankelley.github.io/whalestrike/reference/knot2mps.md)

## Author

Dan Kelley

## Examples

``` r
library(whalestrike)
mps <- seq(0, 10)
plot(mps, mps2knot(mps), xlab = "Speed [m/s]", ylab = "Speed [knots]", type = "l")

```
