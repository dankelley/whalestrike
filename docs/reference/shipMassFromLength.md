# Compute ship displacement from vessel type and length

This is done using formulae in Table 3 of Mayette and Brillant (2026).

## Usage

``` r
shipMassFromLength(type, L)
```

## Arguments

- type:

  a string identifying the ship type, from the list given in ‘Details’.

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

- Alexandra Mayette, Sean W. Brillant. "A regression-based method to
  estimate vessel mass for use in whale-ship strike risk models." PLoS
  One 21(2) e0339760.
  [doi:10.1371/journal.pone.0339760](https://doi.org/10.1371/journal.pone.0339760)

## Author

Dan Kelley, with help from Alexandra Mayette

## Examples

``` r
library(whalestrike)
shipMassFromLength("Tug", 50)/1e3 # 1920.648
#> [1] 1920.648
```
