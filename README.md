# The whalestrike package

<!-- badges: start -->


[![GitHub last commit](https://img.shields.io/github/last-commit/dankelley/whalestrike)](https://img.shields.io/github/last-commit/dankelley/whalestrike)
[![R-CMD-check](https://github.com/dankelley/whalestrike/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dankelley/whalestrike/actions/workflows/R-CMD-check.yaml)
[![status](https://joss.theoj.org/papers/570201320eb0182aa487026819021c50/status.svg)](https://joss.theoj.org/papers/570201320eb0182aa487026819021c50)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)


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
app2()
```

in an R session.  A window will open in your browser.  On the left-hand side is
a region with controllers, many hidden in sub-panels that can be opened by
clicking on V-shaped icons.  Before opening those sub-panels, try adjusting the
ship speed, to see what happens to the plots, especially the one labelled
`Lethality Index`. Think about whether the tendency of variations to the index
are in line with your intuition. Next, open the `Ship` sub-panel, to explore
the result of altering the ship mass.  This sub-panel also has controllers
specifying the geometry of the impact region, and you ought to explore them,
also. Continuing to explore the app's controllers ought to give you a good
indication of what the tool provides. To learn more, try consulting the app's
documentation, provided (a) as the response to typing `? app2` in the R
console, (b) as information in a dialog box that opens when the `Help` button
is clicked, and (c) in a [youtube video](https://youtu.be/kTMl3nXa5A4).

Users who want more control, and who want to deal with the results in numerical
as opposed to graphical form, should try clicking the `Code` button in
`app2()`, to see the R code that runs the simulation outside the app.  The next
step will be to explore the functions used in that code.  Documentation exists
for each of these functions e.g. typing
```R
library(whalestrike)
?strike
```
in an R session will provide information on `strike()`, which is a key function
in the package.

## References

Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. “Assessing the Lethality of Ship
Strikes on Whales Using Simple Biophysical Models.” Marine Mammal
Science, October 12, 2020, mms.12745. ([journal site](https://doi.org/10.1111/mms.12745); 
[ResearchGate Preprint](https://www.researchgate.net/profile/Dan-Kelley-3/publication/344748816_Assessing_the_lethality_of_ship_strikes_on_whales_using_simple_biophysical_models/links/6400e1410cf1030a56678284/Assessing-the-lethality-of-ship-strikes-on-whales-using-simple-biophysical-models.pdf?origin=publicationSearch&_rtd=e30%3D&_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6ImhvbWUiLCJwYWdlIjoic2VhcmNoIiwicG9zaXRpb24iOiJwYWdlSGVhZGVyIn19).)
