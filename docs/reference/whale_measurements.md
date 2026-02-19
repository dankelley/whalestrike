# Get values for various whale measurements

This uses a data frame containing information about several whale
species, compiled by Alexandra Mayette and provided to Dan Kelley as a
personal communication on 2026-01-28.

## Usage

``` r
whale_measurements(species = NULL)
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
