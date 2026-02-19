# Create a function for stress in laminated layers

Denoting unforced layer thickness in the \\i\\ layer as \\l_i\\ and
strain there as \\\epsilon_i=\Delta l_i/l_i\\, we may write the
stress-strain relationship as \$\$\sigma =
a_i\*(exp(b_i\*\epsilon_i)-1)\$\$ for each layer, where it is assumed
that stress \\\sigma\\ is equal across layers. Inverting this yields
\$\$\epsilon_i= ln(1 + \sigma/a_i)/b_i\$\$ where \\ln\\ is the natural
logarithm. Therefore, the change \\\Delta L\\ in the total thickness
\\L=\sum l_i\\ may be written \$\$0 = \Delta L - \sum((l_i/b_i)
ln(1+\sigma/a_i))\$\$. Note that zero-thickness layers are removed from
the calculation, to avoid spurious forces.

## Usage

``` r
stressFromStrainFunction(l, a, b, N = 1000)
```

## Arguments

- l:

  vector of layer thicknesses

- a:

  vector of multipliers

- b:

  vector of e-fold parameters

- N:

  integer specifying how many segments to use in the spline

## Value

A piecewise-linear function, created with
[`approxfun()`](https://rdrr.io/r/stats/approxfun.html), that returns
stress as a function of total strain of the system of compressing
layers. For the purposes of the whale-strike analysis, the strain should
be between 0 and 1, i.e. there is no notion of compressing blubber, etc.
to negative thickness.

## Details

This expression is not easily inverted to get \\\sigma\\ in terms of
\\\Delta L\\ but it may be solved easily for particular numerical
values, using [`uniroot()`](https://rdrr.io/r/stats/uniroot.html).

This is done for a sequence of `N` values of strain \\\epsilon\\ that
range from 0 to 1. Then
[`approxfun()`](https://rdrr.io/r/stats/approxfun.html) is used to
create a piecewise-linear representation of the relationship between
\\\sigma\\ and \\\Delta L\\, which becomes the return value of the
present function. (The purpose of using a piecewise-linear
representation to reduce computation time.)

## See also

Other functions relating to whale characteristics:
[`whaleCompressionForce()`](https://dankelley.github.io/whalestrike/reference/whaleCompressionForce.md),
[`whaleLengthFromMass()`](https://dankelley.github.io/whalestrike/reference/whaleLengthFromMass.md),
[`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md),
[`whaleSkinForce()`](https://dankelley.github.io/whalestrike/reference/whaleSkinForce.md),
[`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md)

Other functions relating to forces:
[`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md),
[`whaleCompressionForce()`](https://dankelley.github.io/whalestrike/reference/whaleCompressionForce.md),
[`whaleSkinForce()`](https://dankelley.github.io/whalestrike/reference/whaleSkinForce.md),
[`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md)

## Author

Dan Kelley

## Examples

``` r
library(whalestrike)
# Set blubber parameters for each layer, to see if
# we recover the raymond2007 data.
param <- parameters(a = rep(1.64e5, 4), b = rep(2.47, 4))
x <- seq(0, 0.5, length.out = 100)
y <- param$stressFromStrain(x)
plot(x, y, type = "l", lwd = 4, col = "gray")
data("raymond2007")
points(raymond2007$strain, raymond2007$stress, col = 2)

```
