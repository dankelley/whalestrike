# Compute whale length from mass

This works by inverting
[`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md)
using [`uniroot()`](https://rdrr.io/r/stats/uniroot.html).

## Usage

``` r
whaleLengthFromMass(M, species = "N. Atl. Right Whale", model = "fortune2012")
```

## Arguments

- M:

  Whale mass (kg).

- species:

  A string indicating the whale species (see
  [`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md)
  for details).

- model:

  Character string specifying the model (see
  [`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md)
  for details).

## Value

Whale length (m).

## References

See
[`whalestrike()`](https://dankelley.github.io/whalestrike/reference/whalestrike.md)
for a list of references.

## See also

[`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md)
is the reverse of this.

Other functions relating to whale characteristics:
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md),
[`whaleCompressionForce()`](https://dankelley.github.io/whalestrike/reference/whaleCompressionForce.md),
[`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md),
[`whaleSkinForce()`](https://dankelley.github.io/whalestrike/reference/whaleSkinForce.md),
[`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md)

## Author

Dan Kelley
