# Compute ship wetted area from mass

Estimate the wetted area of a Cape Islander boat, given the vessel mass.

## Usage

``` r
shipAreaFromMass(ms)
```

## Arguments

- ms:

  Ship mass (kg).

## Value

Estimated area (m^2).

## Details

The method is based on scaling up the results for a single Cape Islander
ship, of displacement 20.46 tonnes, length 11.73m, beam 4.63m, and draft
1.58m, on the assumption that the wetted area is proportional to
\\length\*(2\*draft+beam)\\. This reference area is scaled to the
specified mass, `ms`, by multiplying by the 2/3 power of the mass ratio.

Note that this is a crude calculation meant as a stop-gap measure, for
estimates values of the `Ss` argument to
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md).
It should not be used in preference to inferences made from
architectural drawings of a given ship under study.

## See also

Other functions relating to ship characteristics:
[`shipLength()`](https://dankelley.github.io/whalestrike/reference/shipLength.md),
[`shipMassFromLength()`](https://dankelley.github.io/whalestrike/reference/shipMassFromLength.md),
[`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md)

## Author

Dan Kelley
