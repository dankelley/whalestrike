# Set parameters for a whale strike simulation

Assembles control parameters into a list suitable for passing to
[`strike()`](https://dankelley.github.io/whalestrike/reference/strike.md)
and the functions that it calls. If `file` is provided, then all the
other arguments are read from that source. Note that
[`updateParameters()`](https://dankelley.github.io/whalestrike/reference/updateParameters.md)
may be used to modify the results of `parameters`, e.g. for use in
sensitivity tests.

## Usage

``` r
parameters(
  ms = 45000,
  Ss = NULL,
  Ly = 1.15,
  Lz = 1.15,
  species = "N. Atl. Right Whale",
  lw = 13.7,
  mw = NULL,
  Sw = NULL,
  l = NULL,
  a = NULL,
  b = NULL,
  s = NULL,
  theta = 55,
  Cs = 0.01,
  Cw = 0.0025,
  logistic = list(logStressCenter = 5.38, logStressWidth = 0.349, tau25 = 1e+05, tau50 =
    241000, tau75 = 581000),
  file = NULL
)
```

## Arguments

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

- file:

  Optional name a comma-separated file that holds all of the previous
  values, except `Cs` and `Cw`. If provided, then other parameters
  except `Cs` and `Cw` are ignored, because values are sought from the
  file. The purpose of this is in shiny apps that want to save a
  simulation framework. The file should be saved
  [`write.csv()`](https://rdrr.io/r/utils/write.table.html) with
  `row.names=FALSE`.

## Value

A named list holding the parameters, with defaults and alternatives
reconciled according to the system described above, along with some
items used internally, including `lsum`, which is the sum of the values
in `l`, and `stressFromStrain()`, a function created by
[`stressFromStrainFunction()`](https://dankelley.github.io/whalestrike/reference/stressFromStrainFunction.md)
that computes compression force from engineering strain.

## References

Campbell-Malone, Regina. "Biomechanics of North Atlantic Right Whale
Bone: Mandibular Fracture as a Fatal Endpoint for Blunt Vessel-Whale
Collision Modeling." PhD Thesis, Massachusetts Institute of Technology
and Woods Hole Oceanographic Institution, 2007.
[doi:10.1575/1912/1817](https://doi.org/10.1575/1912/1817) .

Daoust, Pierre-Yves, Emilie L. Couture, Tonya Wimmer, and Laura Bourque.
"Incident Report. North Atlantic Right Whale Mortality Event in the Gulf
of St. Lawrence, 2017." Canadian Wildlife Health Cooperative, Marine
Animal Response Society, and Fisheries and Oceans Canada, 2018.
<https://publications.gc.ca/site/eng/9.850838/publication.html>.

Grear, Molly E., Michael R. Motley, Stephanie B. Crofts, Amanda E. Witt,
Adam P. Summers, and Petra Ditsche. "Mechanical Properties of Harbor
Seal Skin and Blubber - a Test of Anisotropy." Zoology 126 (2018):
137-44.
[doi:10.1016/j.zool.2017.11.002](https://doi.org/10.1016/j.zool.2017.11.002)
.

Raymond, J. J. "Development of a Numerical Model to Predict Impact
Forces on a North Atlantic Right Whale during Collision with a Vessel."
University of New Hampshire, 2007.
<https://scholars.unh.edu/thesis/309/>.

## Author

Dan Kelley

## Examples

``` r
parms <- parameters()
epsilon <- seq(0, 1, length.out = 100) # strain
sigma <- parms$stressFromStrain(epsilon) # stress
plot(epsilon, log10(sigma), xlab = "Strain", ylab = "log10(Stress [MPa])", type = "l")
mtext("Note sudden increase in stress, when bone compression starts")

```
