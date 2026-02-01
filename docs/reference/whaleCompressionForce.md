# Whale compression force

Calculate the total compression stress and force, along with the
thicknesses of skin, blubber, sublayer, and bone. The stress is computed
with the
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md)
function that is created by
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md)
and stored in `para`. the force is computed by multiplying stess by area
computed as the product of `parms$Ly` and `parms$Lz`. Any negative layer
thicknesses are set to zero, as a way to avoid problems with aphysical
engineering compression strains that exceed 1.

## Usage

``` r
whaleCompressionForce(xs, xw, parms)
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

A list containing: `force` (N), the compression-resisting force;
`stress` (Pa), the ratio of that force to the impact area; `strain`, the
total strain, and `compressed`, a four-column matrix (m) with first
column for skin compression, second for blubber compression, third for
sublayer compression, and fourth for bone compression.

## References

See
[`whalestrike()`](https://dankelley.github.io/whalestrike/reference/whalestrike.md)
for a list of references.

## See also

Other functions relating to whale characteristics:
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md),
[`whaleLengthFromMass()`](https://dankelley.github.io/whalestrike/reference/whaleLengthFromMass.md),
[`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md),
[`whaleSkinForce()`](https://dankelley.github.io/whalestrike/reference/whaleSkinForce.md),
[`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md)

Other functions relating to forces:
[`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md),
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md),
[`whaleSkinForce()`](https://dankelley.github.io/whalestrike/reference/whaleSkinForce.md),
[`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md)

## Author

Dan Kelley
