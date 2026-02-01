# Ship displacement in kg based on vessel type and length

This is done using formulae in Table 3 of Mayette and Brillant (2026).

## Usage

``` r
shipMassFromLength(type = NULL, L)
```

## Arguments

- type:

  either (1) a string identifying the ship type, in which case the
  average overall length of the named vessel is returned, or (2) NULL,
  in which case a vector of permitted values of `type` is returned.

- L:

  vessel length in metres.

## Value

`shipMassFromLength` returns ship displacement mass (in kg), according
to Mayette and Brillant (2026) Table 3.

## Details

The formulae used are as follows.

|                       |                  |
|-----------------------|------------------|
| `type`                | Formula          |
| "Bulk Carrier"        | 5.64 \* L^3.06   |
| "Container Ship"      | 86.40 \* L^2.46  |
| "Cruise"              | 97.51 \* L^2.28  |
| "Ferry"               | 25.15 \* L^2.62  |
| "Fishing"             | 0.71 \* L^3.79   |
| "Government/Research" | 2.95 \* L^3.22   |
| "Other"               | 2.64 \* L^3.35   |
| "Passenger"           | 4.32 \* L^3.08   |
| "Pleasure Craft"      | 34.47 \* L^2.68  |
| "Sailing"             | 1.23 \* L^3.53   |
| "Tanker"              | 7.25 \* L^3.03   |
| "Tug"                 | 104.48 \* L^2.51 |

## References

Mayette, Alexandra, and Sean W. Brillant. "A Regression-Based Method to
Estimate Vessel Mass for Use in Whale-Ship Strike Risk Models." PloS One
21, no. 1 (2026): e0339760.
https://doi.org/10.1371/journal.pone.0339760.

## See also

Other functions relating to ship characteristics:
[`shipAreaFromMass()`](https://dankelley.github.io/whalestrike/reference/shipAreaFromMass.md),
[`shipLength()`](https://dankelley.github.io/whalestrike/reference/shipLength.md),
[`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md)

## Author

Dan Kelley, with help from Alexandra Mayette

## Examples

``` r
library(whalestrike)
shipMassFromLength("Tug", 50) / 1e3 # 1920.648
#> [1] 1920.648
```
