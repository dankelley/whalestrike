# Nominal ship length in m

This is based on the "Average LOA in m" column in Table 1 of Mayette and
Brillant (2026).

## Usage

``` r
shipLength(type = NULL)
```

## Arguments

- type:

  either (1) a string identifying the ship type, in which case the
  average overall length of the named vessel is returned, or (2) NULL,
  in which case a data frame containing type and length is returned.

## Value

`shipLength` returns ship length in m, as defined in Mayette and
Brillant (2026).

## References

Mayette, Alexandra, and Sean W. Brillant. "A Regression-Based Method to
Estimate Vessel Mass for Use in Whale-Ship Strike Risk Models." PloS One
21, no. 1 (2026): e0339760.
https://doi.org/10.1371/journal.pone.0339760.

## See also

Other functions relating to ship characteristics:
[`shipAreaFromMass()`](https://dankelley.github.io/whalestrike/reference/shipAreaFromMass.md),
[`shipMassFromLength()`](https://dankelley.github.io/whalestrike/reference/shipMassFromLength.md),
[`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md)

## Author

Dan Kelley, with help from Alexandra Mayette

## Examples

``` r
library(whalestrike)
# An individual length
shipLength("Fishing")
#> [1] 25.9
# A table of lengths
shipLength()
#>                   type length
#> 1         Bulk Carrier  221.0
#> 2       Container Ship  240.7
#> 3               Cruise  219.1
#> 4                Ferry  100.1
#> 5              Fishing   25.9
#> 6  Government/Research   63.7
#> 7                Other   83.3
#> 8            Passenger   23.2
#> 9       Pleasure Craft   42.2
#> 10             Sailing   20.3
#> 11              Tanker  198.4
#> 12                 Tug   44.1
```
