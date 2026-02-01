# Whale Drag Force

Compute the retarding force of water on the whale, based on a drag law
\\(1/2)\*rho\*Cw\*A\*vw^2\\ where `rho` is 1024 (kg/m^3), `Cw` is
`parms$Cw` and `A` is `parms$Sw`.

## Usage

``` r
whaleWaterForce(vw, parms)
```

## Arguments

- vw:

  Whale velocity (m/s).

- parms:

  A named list holding model parameters, created by
  [`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md).

## Value

Water drag force (N).

## See also

Other functions relating to whale characteristics:
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md),
[`whaleCompressionForce()`](https://dankelley.github.io/whalestrike/reference/whaleCompressionForce.md),
[`whaleLengthFromMass()`](https://dankelley.github.io/whalestrike/reference/whaleLengthFromMass.md),
[`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md),
[`whaleSkinForce()`](https://dankelley.github.io/whalestrike/reference/whaleSkinForce.md)

Other functions relating to forces:
[`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md),
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md),
[`whaleCompressionForce()`](https://dankelley.github.io/whalestrike/reference/whaleCompressionForce.md),
[`whaleSkinForce()`](https://dankelley.github.io/whalestrike/reference/whaleSkinForce.md)

## Author

Dan Kelley
