# Skin force

The ship-whale separation is used to calculate the deformation of the
skin. The parameters of the calculation are `parms$Ly` (impact area
width, m), `parms$Lz` (impact area height, in m), `parms$Ealpha` (skin
elastic modulus in Pa), `parms$alpha` (skin thickness in m), and
`parms$theta` (skin bevel angle degrees, measured from a vector normal
to undisturbed skin).

## Usage

``` r
whaleSkinForce(xs, xw, parms)
```

## Arguments

- xs:

  Ship position (m).

- xw:

  Whale position (m).

- parms:

  A named list holding model parameters, created by
  [`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md).

## Value

A list containing `force`, the normal force (N), along with `sigmay` and
`sigmaz`, which are stresses (Pa) in the y (beam) and z (draft)
directions.

## References

See
[`whalestrike()`](https://dankelley.github.io/whalestrike/reference/whalestrike.md)
for a list of references.

## See also

Other functions relating to whale characteristics:
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md),
[`whaleCompressionForce()`](https://dankelley.github.io/whalestrike/reference/whaleCompressionForce.md),
[`whaleLengthFromMass()`](https://dankelley.github.io/whalestrike/reference/whaleLengthFromMass.md),
[`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md),
[`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md)

Other functions relating to forces:
[`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md),
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md),
[`whaleCompressionForce()`](https://dankelley.github.io/whalestrike/reference/whaleCompressionForce.md),
[`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md)

## Author

Dan Kelley
