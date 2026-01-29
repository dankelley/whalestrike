# Whale measurements

`whale_measurements` is a data frame containing information about
several whale species, compiled by Alexandra Mayette and provided to Dan
Kelley as a personal communication on 2026-01-28.

## Details

The columns are follows.

- `name` species name, as used in this package.

- `Species` proper species name, not used in this package.

- `length` whale length in metres.

- `girth` whale girth in centimetres.

- `bone` whale bone thickess in centimeters.

- `sublayer` thickness of sublayer in centimeters; this was called
  muscle in the Mayette document

- `blubber` whale blubber thickness in centimeters.

- `skin` whale skin thicness in centimetres.

## Examples

``` r
library(whalestrike)
data(whale_measurements)
whale_measurements
#>            name                    Species length girth bone muscle blubber
#> 1          Blue      Balaenoptera musculus   21.9  1219 8.65  168.7     7.6
#> 2           Fin      Balaenoptera physalus   16.9   546 7.20   71.4     5.8
#> 3          Gray      Eschrichtius robustus   12.5   685 6.80   84.7     9.7
#> 4      Humpback      Megaptera novaengliae   11.6  1008 7.10  137.2     8.2
#> 5         Minke Balaenoptera acutorostrata    6.7   268 3.90   31.3     3.3
#> 6 N. Atl. Right        Eubalaena glacialis   13.8  1030 7.15  132.5    16.3
#> 7           Sei      Balaenoptera borealis   13.4   474 4.65   63.1     4.4
#> 8         Sperm     Physeter macrocephalus   11.8   686 8.20   81.0    11.5
#>   skin
#> 1  0.4
#> 2  0.5
#> 3  1.0
#> 4  0.9
#> 5  0.3
#> 6  0.9
#> 7  0.2
#> 8  0.4
```
