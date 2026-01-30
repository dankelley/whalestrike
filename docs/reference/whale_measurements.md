# Whale measurements

`whale_measurements` is a data frame containing information about
several whale species, compiled by Alexandra Mayette and provided to Dan
Kelley as a personal communication on 2026-01-28.

## Details

There are two additions to `whale_measurements` that are not in
Mayette's table. These are `"Pac. Right"` and `"Bryde"`. For these,
Mayette has suggesting using values for the `"N. Atl. Right"` and
`"Sei"` cases, respectively, as conditional estimates for use in this
package.

The columns are follows.

- `name` species name, as used in e.g. whaleMassFromLength().

- `Species` proper species name. (This is not used in this package.)

- `length` whale length in metres.

- `bone` whale bone thickness (measured from the centre to the
  sublayer).

- `sublayer` thickness of sublayer in centimeters; this was called
  muscle in the Mayette document

- `blubber` whale blubber thickness in centimeters.

- `skin` whale skin thickness in centimetres.

## Examples

``` r
library(whalestrike)
data(whale_measurements)
whale_measurements
#>             name                    Species length bone muscle blubber skin
#> 1           Blue      Balaenoptera musculus   21.9 17.3  168.7     7.6  0.4
#> 2          Bryde        Balaenoptera brydei   13.4  9.3   63.1     4.4  0.2
#> 3            Fin      Balaenoptera physalus   16.9 14.4   71.4     5.8  0.5
#> 4           Gray      Eschrichtius robustus   12.5 13.6   84.7     9.7  1.0
#> 5       Humpback      Megaptera novaengliae   11.6 14.2  137.2     8.2  0.9
#> 6          Minke Balaenoptera acutorostrata    6.7  7.8   31.3     3.3  0.3
#> 7  N. Atl. Right        Eubalaena glacialis   13.8 14.3  132.5    16.3  0.9
#> 8     Pac. Right         Eubalaena japonica   13.8 14.3  132.5    16.3  0.9
#> 9            Sei      Balaenoptera borealis   13.4  9.3   63.1     4.4  0.2
#> 10         Sperm     Physeter macrocephalus   11.8 16.4   81.0    11.5  0.4
```
