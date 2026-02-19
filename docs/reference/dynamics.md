# Dynamical law

This function handles Newton's second law, which is the dynamical law
that relates the accelerations of whale and ship to the forces upon
each. It is used by
[`strike()`](https://dankelley.github.io/whalestrike/reference/strike.md),
as the latter integrates the acceleration equations to step forward in
time through the simulation of a whale-strike event. Thus, `dynamics()`
is a core function of this package. The code is very simple, because the
forces are determined by other functions, as described in the “Details”
section.

## Usage

``` r
dynamics(t, y, parms)
```

## Arguments

- t:

  time (s).

- y:

  model state, a vector containing ship position `xs` (m), ship speed
  `vs` (m/s), whale position `xw` (m), and whale speed `vw` (m/s).

- parms:

  A named list holding model parameters, created by
  [`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md).

## Details

Given a present state (defined by the positions and velocities of ship
and whale) at the present time, apply Newton's second law to find the
time derivatives of that state. Forces are determined with
[`whaleCompressionForce()`](https://dankelley.github.io/whalestrike/reference/whaleCompressionForce.md),
[`whaleSkinForce()`](https://dankelley.github.io/whalestrike/reference/whaleSkinForce.md),
[`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md),
[`whaleWaterForce()`](https://dankelley.github.io/whalestrike/reference/whaleWaterForce.md),
while engine force (assumed constant over the course of a collision) is
computed from initial
[`shipWaterForce()`](https://dankelley.github.io/whalestrike/reference/shipWaterForce.md).
Whale and ship masses are set by
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md),
which also sets up areas, drag coefficients, etc.

## References

See
[`whalestrike()`](https://dankelley.github.io/whalestrike/reference/whalestrike.md)
for a list of references.

## Author

Dan Kelley
