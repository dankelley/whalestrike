# Find maximum Lethality Index during a strike

This works by finding the maximum Lethality Index encountered during a
simulation created by calling
[`strike()`](https://dankelley.github.io/whalestrike/reference/strike.md),
and so it is important to use a detailed setting for the output times.
In the example, the results are reported every 0.7/200 seconds (i.e. 3.5
milliseconds), which is likely sufficient (see the example, where a plot
is used for this assessment).

## Usage

``` r
maximumLethalityIndex(strike)
```

## Arguments

- strike:

  the value returned by a call to strike.

## See also

Other functions dealing with Whale Lethality index:
[`lethalityIndexFromStress()`](https://dankelley.github.io/whalestrike/reference/lethalityIndexFromStress.md),
[`stressFromLethalityIndex()`](https://dankelley.github.io/whalestrike/reference/stressFromLethalityIndex.md)

## Author

Dan Kelley, wrapping code provided by Alexandra Mayette

## Examples

``` r
library(whalestrike)
t <- seq(0, 0.7, length.out = 200)
state <- list(xs = -2, vs = knot2mps(10), xw = 0, vw = 0)
parms <- parameters()
s <- strike(t, state, parms)
# Compute the desired value and (for context) show it on a plot
maximumLethalityIndex(s)
#> [1] 0.6931811
# For context, this is how this can be done "by hand"
max(lethalityIndexFromStress(s[["WCF"]][["stress"]]))
#> [1] 0.6931811
# Show the maximum on a plot (see also the plot title)
plot(s, which = "lethality index")
abline(h=maximumLethalityIndex(s), col=2)

```
