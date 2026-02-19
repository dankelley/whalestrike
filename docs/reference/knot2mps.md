# Convert a speed in knots to a speed in m/s

See also
[`mps2knot()`](https://dankelley.github.io/whalestrike/reference/mps2knot.md),
which is the inverse of this function.

## Usage

``` r
knot2mps(knot)
```

## Arguments

- knot:

  Speed in knots.

## Value

Speed in m/s.

## See also

Other functions dealing with units:
[`mps2knot()`](https://dankelley.github.io/whalestrike/reference/mps2knot.md)

## Author

Dan Kelley

## Examples

``` r
library(whalestrike)
knots <- seq(0, 20)
plot(knots, knot2mps(knots), xlab = "Speed [knots]", ylab = "Speed [m/s]", type = "l")

```
