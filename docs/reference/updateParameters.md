# Update parameters

`updateParameters()` is used to alter one or more components of an
existing object of type `"parameters"` that was created by
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md).
This can be useful for e.g. sensitivity tests (see “Details”).

## Usage

``` r
updateParameters(
  original,
  ms,
  Ss,
  Ly,
  Lz,
  species,
  lw,
  mw,
  Sw,
  l,
  a,
  b,
  s,
  theta,
  Cs,
  Cw,
  logistic,
  debug = 0
)
```

## Arguments

- original:

  An object of class `"parameters"`, as created by
  [`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md)
  and perhaps later altered by previous calls to `updateParameters()`.

- ms:

  Ship mass (kg).

- Ss:

  Ship wetted area (m^2). This, together with `Cs`, is used by
  [`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md)
  to estimate ship drag force. If `Ss` is not given, then an estimate is
  made by calling
  [`shipAreaFromMass()`](https://dankelley.github.io/whalestrike/reference/shipAreaFromMass.md)
  with the provided value of `ms`.

- Ly:

  Ship impact horizontal extent (m); defaults to 1.15m if not specified,
  based on an analysis of the shape of the bow of typical coastal
  fishing boats of the Cape Islander variety.

- Lz:

  Ship impact vertical extent (m); defaults to 1.15m if not specified,
  based on the same analysis as for Ly.

- species:

  a string indicating the whale species. For the permitted values, see
  [`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md).
  (The `species` value can also set the `lw` and `l` values, as noted in
  their portions of this documentation.)

- lw:

  either (1) whale length in metres or (2) the string `"from_species"`.
  If the latter, then the length is determined from
  [`whaleMeasurements()`](https://dankelley.github.io/whalestrike/reference/whaleMeasurements.md).
  In either case, the length is used by
  [`whaleAreaFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleAreaFromLength.md)
  to calculate area, which is needed for the water drag calculation done
  by
  [`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md).

- mw:

  either (1) the whale mass in kg or (2) NULL. In the latter case, the
  mass is calculated from whale length, using
  [`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md)
  with `type="wetted"`.

- Sw:

  either (1) the whale surface area in m^2 or (2) NULL. If the latter
  case, the area is calculated from whale length using
  [`whaleAreaFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleAreaFromLength.md).

- l:

  either (1) a numerical vector of length 4 that indicates the
  thicknesses in metres of skin, blubber, sublayer and bone; (2) NULL to
  set these four values to 0.025, 0.16, 1.12, and 0.1; or (3) the string
  `"from_species"`, in which case these four values are determined by
  calling
  [`whaleMeasurements()`](https://dankelley.github.io/whalestrike/reference/whaleMeasurements.md).
  The default skin thickness of 0.025 m represents the 0.9-1.0 inch
  value stated in Section 2.2.3 of Raymond (2007). The blubber default
  of 0.16 m is a rounded average of the values inferred by whale
  necropsy, reported in Appendix 2 of Daoust et al., 2018. The sublayer
  default of 1.12 m may be reasonable at some spots on the whale body.
  The bone default of 0.1 m may be reasonable at some spots on the whale
  body. The sum of these default values, 1.40 m, is a whale radius that
  is consistent with a half-circumference of 4.4 m, reported in Table
  2.2 of Raymond (2007). Note, however, that these values are not
  identical to those found in `whaleMeasurements`.

- a, b:

  Numerical vectors of length 4, giving values to use in the
  stress-strain law `stress=a*(exp(b*strain)-1)`, where `a` is in Pa and
  `b` is unitless. By construction, `a*b` is the local modulus at low
  strain (i.e. at low `b*strain` values), and that `b` is the efolding
  scale for nonlinear increase in stress with strain. This exponential
  relationship has been mapped out for whale blubber, using a curve fit
  to Figure 2.13 of Raymond (2007), and these values are used for the
  second layer (blubber); see the documentation for the
  [raymond2007](https://dankelley.github.io/whalestrike/reference/raymond2007.md)
  dataset, to see for how that fit was done. If not provided, `a`
  defaults to `c(17.8e6/0.1, 1.58e5, 1.58e5, 8.54e8/0.1)` and `b`
  defaults to `c(0.1, 2.54, 2.54, 0.1)`. The skin defaults are set up to
  give a linear shape (since `b` is small) with the `a*b` product being
  17.8e6 Pa, which is the adult-seal value given in Table 3 of Grear et
  al. (2017). The blubber defaults are from a regression of the
  stress-strain relationship shown in Figure 2.13 of Raymond (2007). The
  sublayer defaults are set to match those of blubber, lacking any other
  information. The bone default for `b` is small, to set up a linear
  function, and `a*b` is set to equal 8.54e8 Pa, given in Table 2.3 of
  Raymond (2007) and Table 4.5 of Campbell-Malone (2007).

- s:

  Numerical vector of length 4, giving the ultimate strengths (Pa) of
  skin, blubber, sublayer, and bone, respectively. If not provided, the
  value is set to `1e6 * c(19.600,0.255,0.255,22.900)` with reasoning as
  follows. The skin default of 19.6 MPa is a rounded value from Table 3
  of Grear et al. (2018) for adult seal skin strength at an orientation
  of 0 degrees. The blubber and sublayer values were chosen as the
  central point of a logistic fit of whale collision damage to maximal
  stress during a default impact simulation. (For comparison, a strength
  of 0.437 MPa may be inferred by multiplying Raymond's (2007) Figure
  2.13 elastic modulus of 0.636 MPa by the ratio 0.97/1.41 determined
  for adult seal strength/modulus, as reported in Table 3 of Grear et
  al. (2018).) The bone default o 22.9 MPa is from Table 2.3 of
  Raymond (2007) and Table 4.5 of Campbell-Malone (2007).

- theta:

  Whale skin deformation angle (deg); defaults to 55 degrees, if not
  supplied, because that angle produces a good match to Raymond's (2007)
  Figure 6.1 for the total force as a function of vessel speed, for
  large vessels. Note that the match works almost as well in the range
  50 deg to 70 deg.

- Cs:

  Drag coefficient for ship (dimensionless), used by
  [`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md)
  to estimate ship drag force. Defaults to 1e-2, which is 4 times the
  frictional coefficient of 2.5e-3 inferred from Figure 4 of Manen and
  van Oossanen (1988), assuming a Reynolds number of 5e7, computed from
  speed 5m/s, lengthscale 10m and viscosity 1e-6 m^2/s. The factor of 4
  is under the assumption that frictional drag is about a quarter of
  total drag. The drag force is computed with
  [`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md).

- Cw:

  Drag coefficient for whale (dimensionless), used by
  [`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md)
  to estimate whale drag force. Defaults to 2.5e-3, for Reynolds number
  2e7, computed from speed 2 m/s, lengthscale 5m which is chosen to be
  between radius and length, and viscosity 1e-6 m^2/s. The drag force is
  computed with
  [`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md).

- logistic:

  a [list](https://rdrr.io/r/base/list.html) containing
  `logStressCenter` and `logStressWidth`, which define an empirical
  logistic fit of an index of whale injury in observed strikes (ranging
  from 0 for no injury to 1 for fatal injury), as a function of the
  base-10 logarithm of compressive stress, as well as `tau25`, `tau50`
  and `tau75`, which are the stresses in that fit that yield index
  values of 0.25, 0.50 and 0.75, respectively; these values set colour
  boundaries in
  [`plot.strike()`](https://dankelley.github.io/whalestrike/reference/plot.strike.md)
  plots that have `which="threat"`.

- debug:

  Integer indicating debugging level, 0 for quiet operation and higher
  values for more verbose monitoring of progress through the function.

## Value

A named list holding the items of the same name as those in the list
returned by
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md).

## Details

Two important differences between argument handling in
`updateParameters()` and
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md)
should be kept in mind.

First, `updateParameters()` does not check its arguments for feasible
values. This can lead to bad results when using
[`strike()`](https://dankelley.github.io/whalestrike/reference/strike.md),
which is e.g. expecting four layer thicknesses to be specified, and also
that each thickness is positive.

Second, `updateParameters()` does not perform ancillary actions that
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md)
performs, with regard to certain interlinking argument values. Such
actions are set up for whale length and ship mass, which are
easily-observed quantities from other quantities can be estimated using
`whalestrike` functions. If `lw` (whale length) is supplied to
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md)
without also supplying `mw` (whale mass), then
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md)
uses
[`whaleMassFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleMassFromLength.md)
to infer `mw` from `lw`. The same procedure is used to infer `Sw` if it
is not given, using
[`whaleAreaFromLength()`](https://dankelley.github.io/whalestrike/reference/whaleAreaFromLength.md).
Similarly,
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md)
uses
[`shipAreaFromMass()`](https://dankelley.github.io/whalestrike/reference/shipAreaFromMass.md)
to compute `Ss` (ship area) from `ms` (ship mass), if the `Ss` argument
is not given. Importantly, these three inferences are *not* made by
`updateParameters()`, which alters only those values that are supplied
explicitly. It is easy to supply those values, however; for example,

    parms <- updateParameters(PARMS, lw=1.01 * PARMS$lw))
    parms <- updateParameters(parms, mw=whaleMassFromLength(parms$lw))
    parms <- updateParameters(parms, Sw=whaleAreaFromLength(parms$lw))

modifies a base state stored in `PARMS`, increasing whale length by 1%
and then increasing whale mass and area accordingly. This code block is
excerpted from a sensitivity test of the model, in which

    parms <- updateParameters(PARMS, ms=1.01 * PARMS$ms)
    parms <- updateParameters(parms, Ss=shipAreaFromMass(parms$ms))

was also used to perturb ship mass (and inferred area).

## References

Daoust, Pierre-Yves, Emilie L. Couture, Tonya Wimmer, and Laura Bourque.
"Incident Report. North Atlantic Right Whale Mortality Event in the Gulf
of St. Lawrence, 2017." Canadian Wildlife Health Cooperative, Marine
Animal Response Society, and Fisheries and Oceans Canada, 2018.
<https://publications.gc.ca/site/eng/9.850838/publication.html>.

## Author

Dan Kelley
