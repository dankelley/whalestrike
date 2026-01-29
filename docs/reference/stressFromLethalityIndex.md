# Compute stress, based on lethality index

The model used for this is the logistic model, fitting observed
injury/lethality statistics to the base-10 logarithm of the maximum
compression stress during a simulated impact event.

## Usage

``` r
stressFromLethalityIndex(injury)
```

## Arguments

- injury:

  numerical value or vector, giving threat of injury (in range 0 to 1).

## Value

whale compression stress, in Pascals.

## See also

Other functions dealing with Whale Lethality index:
[`lethalityIndexFromStress()`](https://dankelley.github.io/whalestrike/reference/lethalityIndexFromStress.md),
[`maximumLethalityIndex()`](https://dankelley.github.io/whalestrike/reference/maximumLethalityIndex.md)

## Author

Dan Kelley

## Examples

``` r
stressFromLethalityIndex(0.5) # approx. 254000 Pa, i.e. parameters()$logistic$tau50
#> [1] 239883.3
```
