# Get values for various whale measurements

This uses a data frame containing information about several whale
species, compiled by Alexandra Mayette and provided to Dan Kelley as a
personal communication on 2026-01-28.

## Usage

``` r
whaleMeasurements(species = NULL)
```

## Arguments

- species:

  either (1) the name of a species or (2) NULL. In the first case, the
  table is consulted to find a row with the given species name, and that
  row is returned. In the second case, the whole table is returned.

## Value

The return value contains

- `name` species name, as used in e.g. whaleMassFromLength().

- `Species` proper species name. (This is not used in this package.)

- `length` whale length in metres.

- `bone` whale bone thickness in metres, measured from the centre to the
  sublayer.

- `sublayer` thickness of sublayer in meters; this was called `muscle`
  in Mayette's document.

- `blubber` whale blubber thickness in meters.

- `skin` whale skin thickness in metres.

## Details

There are two species in the table that are not in Mayette's table.
These are `"Pac. Right Whale"` and `"Bryde Whale"`. For these, Mayette
has suggesting using values for the `"N. Atl. Right Whale"` and
`"Sei Whale"` cases, respectively, as conditional estimates for use in
this package.

## Author

Dan Kelley, using data and advice from Alexandra Mayette

## Examples

``` r
library(whalestrike)
# All species in database
whaleMeasurements()
#>                species                 properName length  bone sublayer blubber
#> 1           Blue Whale      Balaenoptera musculus   21.9 0.173    1.687   0.076
#> 2          Bryde Whale        Balaenoptera brydei   13.4 0.093    0.631   0.044
#> 3            Fin Whale      Balaenoptera physalus   16.9 0.144    0.714   0.058
#> 4           Gray Whale      Eschrichtius robustus   12.5 0.136    0.847   0.097
#> 5       Humpback Whale      Megaptera novaengliae   11.6 0.142    1.372   0.082
#> 6          Minke Whale Balaenoptera acutorostrata    6.7 0.078    0.313   0.033
#> 7  N. Atl. Right Whale        Eubalaena glacialis   13.8 0.143    1.325   0.163
#> 8     Pac. Right Whale         Eubalaena japonica   13.8 0.143    1.325   0.163
#> 9            Sei Whale      Balaenoptera borealis   13.4 0.093    0.631   0.044
#> 10         Sperm Whale     Physeter macrocephalus   11.8 0.164    0.810   0.115
#>     skin
#> 1  0.004
#> 2  0.002
#> 3  0.005
#> 4  0.010
#> 5  0.009
#> 6  0.003
#> 7  0.009
#> 8  0.009
#> 9  0.002
#> 10 0.004
# A particular species
whaleMeasurements("N. Atl. Right Whale")
#>               species          properName length  bone sublayer blubber  skin
#> 7 N. Atl. Right Whale Eubalaena glacialis   13.8 0.143    1.325   0.163 0.009
```
