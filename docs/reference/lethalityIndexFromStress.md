# Compute lethality index, based on compression stress

The model used for this is the logistic model, fitting observed
injury/lethality statistics to the base-10 logarithm of the maximum
compression stress during a simulated impact event.

## Usage

``` r
lethalityIndexFromStress(stress)
```

## Arguments

- stress:

  numerical value or vector, giving whale compression stress in Pascals.

## Value

threat of injury (in range 0 to 1)

## See also

Other functions dealing with Whale Lethality index:
[`maximumLethalityIndex()`](https://dankelley.github.io/whalestrike/reference/maximumLethalityIndex.md),
[`stressFromLethalityIndex()`](https://dankelley.github.io/whalestrike/reference/stressFromLethalityIndex.md)

## Author

Dan Kelley

## Examples

``` r
lethalityIndexFromStress(parameters()$logistic$tau50) # approx. 0.5
#> [1] 0.5014449
```
