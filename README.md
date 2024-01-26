# The whalestrike package

<!-- badges: start -->

[![GitHub last commit](https://img.shields.io/github/last-commit/dankelley/whalestrike)](https://img.shields.io/github/last-commit/dankelley/whalestrike)
[![R-CMD-check](https://github.com/dankelley/whalestrike/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dankelley/whalestrike/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Purpose

This package provides tools for simulating the collisions of ships with whales,
using a simplified dynamical structure involving point masses separated by
compressible materials. Along with functions for computations of forces,
deformations, accelerations, etc., the package provides a easy-to-use GUI tool
that makes it easy to set up some common simulation scenarios, and to display
the results in graphical form.

To learn more about the scientific background, and to see the results of some
detailed computations placed in the context of a database of observed strikes,
see Kelley et al. (2020).

## Installation

The package is not yet available on CRAN, and must be installed from source.
This can be done either by downloading the source and building it locally, or
by the simpler method of typing
```R
# install.package("remotes")
remotes::install_github("dankelley/whalestrike", ref="main")
```
in an R session.  (Uncomment the first line, if the `remotes`
package is not yet installed on your machine.)

## Usage

Most users will find that the GUI application is a good way to learn about the
package.  To start this, type
```R
library(whalestrike)
app()
```
in an R session.  A window will open, which shows some information at the top.
Below that several sliders, buttons, and pull-down menus provide ways to control
the simulation and the graphical representation of the results. Use `?app` to
learn more about the process.

Users who want more control, and who want to deal with the results in numerical
as opposed to graphical form, should then start exploring the `strike()`
function.  To learn how it works, type
```R
library(whalestrike)
?strike
```
in an R session.

## References

Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. “Assessing the Lethality of Ship
Strikes on Whales Using Simple Biophysical Models.” Marine Mammal
Science, October 12, 2020, mms.12745. (https://doi.org/10.1111/mms.12745).

