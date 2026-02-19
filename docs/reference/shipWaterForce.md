# Ship Drag Force

Compute the retarding force of water on the ship, based on a drag law
\\(1/2)\*rho\*Cs\*A\*vs^2\\ where `rho` is water density taken to be
1024 (kg/m^3), `Cs` is drag coefficient stored in `parms` and `A` is
area, also stored in `parms, and `vs\` is the ship speed (m/s).

## Usage

``` r
shipWaterForce(vs, parms)
```

## Arguments

- vs:

  Ship speed in m/s. (Consider using
  [`knot2mps()`](https://dankelley.github.io/whalestrike/reference/knot2mps.md)
  if you prefer to think of speeds in knots.)

- parms:

  A named list holding model parameters, created by
  [`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md).

## Value

Water drag force (N).

## See also

Other functions relating to ship characteristics:
[`shipAreaFromMass()`](https://dankelley.github.io/whalestrike/reference/shipAreaFromMass.md),
[`shipLength()`](https://dankelley.github.io/whalestrike/reference/shipLength.md),
[`shipMassFromLength()`](https://dankelley.github.io/whalestrike/reference/shipMassFromLength.md)

Other functions relating to forces:
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md),
[`whaleCompressionForce()`](https://dankelley.github.io/whalestrike/reference/whaleCompressionForce.md),
[`whaleSkinForce()`](https://dankelley.github.io/whalestrike/reference/whaleSkinForce.md),
[`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md)

## Author

Dan Kelley
