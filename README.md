# The whalestrike package

<!-- badges: start -->


[![GitHub last commit](https://img.shields.io/github/last-commit/dankelley/whalestrike)](https://img.shields.io/github/last-commit/dankelley/whalestrike)
[![R-CMD-check](https://github.com/dankelley/whalestrike/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dankelley/whalestrike/actions/workflows/R-CMD-check.yaml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.06473/status.svg)](https://doi.org/10.21105/joss.06473)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11372537.svg)](https://doi.org/10.5281/zenodo.11372537)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)


<!-- badges: end -->

## Purpose

This package provides tools for simulating the collisions of ships
with whales, using a simplified dynamical structure involving point
masses separated by compressible materials. Along with functions for
computations of forces, deformations, accelerations, etc., the package
provides a easy-to-use GUI tool that makes it easy to set up some
common simulation scenarios, and to display the results in graphical
form.

To learn more about the scientific background, and to see the results
of some detailed computations placed in the context of a database of
observed strikes, see Kelley et al. (2021). For more information on
the coding, see Kelley (2024).

## Installation

The package is not yet available on CRAN, and must be installed from
source. This can be done either by downloading the source and building
it locally, or by the simpler method of typing
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

in an R session.  A window will open in your browser.  On the
left-hand side is a region with controllers, many hidden in sub-panels
that can be opened by clicking on V-shaped icons.  Before opening
those sub-panels, try adjusting the ship speed, to see what happens to
the plots, especially the one labelled `Lethality Index`. Think about
whether the tendency of variations to the index are in line with your
intuition. Next, open the `Ship` sub-panel, to explore the result of
altering the ship mass.  This sub-panel also has controllers
specifying the geometry of the impact region, and you ought to explore
them, also. Continuing to explore the app's controllers ought to give
you a good indication of what the tool provides. To learn more, try
consulting the app's documentation, provided (a) as the response to
typing `? app2` in the R console, (b) as information in a dialog box
that opens when the `Help` button is clicked, and (c) in a [youtube
video](https://youtu.be/kTMl3nXa5A4) and a brief [followup youtube
video](https://youtu.be/f8nHGikb9ug) that illustrates an alteration
made after a helpful comment from a review of a Journal of Open-Source
Software manuscript about the package.

Users who want more control, and who want to deal with the results in
numerical as opposed to graphical form, should try clicking the `Code`
button in `app2()`, to see the R code that runs the simulation outside
the app.  The next step will be to explore the functions used in that
code.  Documentation exists for each of these functions e.g. typing
```R
library(whalestrike)
?strike
```
in an R session will provide information on `strike()`, which is a key
function in the package.

## References


* Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. "Assessing
  the Lethality of Ship Strikes on Whales Using Simple Biophysical
  Models." Marine Mammal Science 37, no. 1 (January 2021): 251–67.
  https://doi.org/10.1111/mms.12745. ([journal
  site](https://doi.org/10.1111/mms.12745); [ResearchGate
  Preprint](https://www.researchgate.net/publication/344748816_Assessing_the_lethality_of_ship_strikes_on_whales_using_simple_biophysical_models))

* Kelley, Dan E. "Whalestrike: an R package for simulating ship
  strikes on whales". Submitted to Journal of Open Source Software,
  2024-01-29.
