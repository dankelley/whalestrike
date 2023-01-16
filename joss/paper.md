---
title: 'whalestrike: An R package for simulating ship strikes on whales'
tags:
  - R
  - oceanography
  - whales
  - ship strikes
  - marine policy
authors:
  - name: Dan E. Kelley
    orcid: 0000-0001-7808-5911
    affiliation: 1
affiliations:
 - name: Department of Oceanography, Dalhousie University, Halifax, Nova Scotia, Canada
   index: 1
citation_author: Dan Kelley
date: 1 February 2023
year: 2023
bibliography: paper.bib
---

# Summary

The collision of ships with whales can result in serious injury or death of the
animal.  This is particularly concerning for endangered species such as the
North Atlantic right whale. Speed reduction policies have been developed, but
simple considerations suggest that other factors, such as ship mass and contact
area, should also be taken into account. An R package called `whalestrike` has
been developed to address these issues. It has been used by
@kelley_assessing_2020 in the development of a biomechanically based criterion
for the lethality of ship strikes. However, this is just one purpose to which
the package could be put. Accordingly, the goal of the present paper is to
introduce readers to the `whalestrike` code, as a way to encourage its wider
usage and development, by researchers and perhaps also by policy makers.

# Statement of need

The collision of a ship with a whale can result in serious injury or death of
the animal, posing significant threats for endangered species
[@laist_collisions_2001]. Of particular concern is the North Atlantic right
whale (*Eubalaena glacialis*), a "critically Endangered" species
[@iucn_eubalaena_2020] with a world population estimated to be just $336\pm14$
in 2021 [@pettis_north_2022], down from $483$ in 2010 [@pace_state-space_2017].
Necropsies reveal that ship collisions account for more than half of right
whale deaths [@campbell-malone_gross_2008-1]. Motivated by such studies,
efforts have been made in recent years to mitigate collision risk by imposing
speed restrictions on ships, and evidence of successful results [e.g.
@conn_vessel_2013] has led to marine policy changes, including both static and
dynamic zones of speed restriction [e.g. @transport_canada_protecting_2022].

Still, it seems unwise to measure the success of speed restrictions by counting
dead or injured whales, given the low numbers alive today. In addition, basic
reasoning reveals that speed cannot be the only factor. Variations in ship
mass, prow shape, etc. should also be considered, making for a multifactorial
problem that requires even more data to achieve statistical reliability.

With this in mind, a group of us undertook a study using a numerical model of
the biophysical dynamics of collisions between vessels and whales, relying on
published records of whale injury and death to calibrate a lethality criterion
as a function of vessel mass, prow geometry, etc., in addition to speed
[@kelley_assessing_2020]. The model is expressed in the R language, which is
familiar to many marine biologists and which provides a wide scope of
statistical tools that might be employed for followup work. The purpose of the
present paper is not to recapitulate the results of @kelley_assessing_2020,
but rather to introduce readers to the model code, in the hopes that they might
suggest extensions or use it for new applications.

# Model formulation

A desire to produce a GUI-based tool permitting easy exploration of various
model scenarios led to a decision to create a simplified model that yields
results quickly, as opposed to a much more computationally-expensive finite
element model that might account more accurately for the deformation of whale
flesh [e.g. that of @raymond_development_2007].  The model ignored ship
deformation upon impact, and considered whale deformation to occur only in a
specified impact area dictated by the shape of the ship's prow.  A layered
model was used for the whale, with skin covering blubber, which in turn covered
what we called a sub-layer (representing a combination of muscle and organs),
and with bone at the core. Thicknesses and material properties for each layer
were taken from the literature, and adjusting these properties permits
simulation of strikes on different species, or at different body locations.
Forces associated with the skin deformation include both extension forces and
compression forces, while the only compression was considered for the other
layers. Lacking reliable information on failure limits for these biomaterials,
critical values for stresses were inferred by reference to published results of
the damage experienced in documented ship strikes. Readers seeking more
information on the parameterizations employed in the model are directed to
@kelley_assessing_2020.

The model mechanics are simple, with acceleration and forces being linked via
Newton's second law. Numerical integration of acceleration is done with the
`lsoda()` function of the `deTools` package.  A first integration yields
velocities, which are used in computing water drag.  A second integration
yields the relative position ship and whale, from which extension and
compression forces can be computed using a nonlinear stress-strain relationship
along with the contact area.  Aspects of each dynamical element are distilled
into the nearly 30 parameters of the model.  Although a great deal of effort
has been put into formulating these parameters appropriately (mainly for
application to right whales), `whalestrike` offers a simple way for users to
adjust each of them, if needed for new applications.

# Package installation and use

Being developed on github.com/dankelley/whalestrike, the package may be installed by typing
```R
library(remotes)
install_github("dankelley/whalestrike", ref="main")
```
in an R console^[There are plans to submit `whalestrike` to the Comprehensive
R Archive Network [@cran_comprehensive_2023].].  Once this is done, the user
has access to over 20 functions and their documentation.  The latter makes
frequent reference to the scientific literature, since the package is aimed at
scientists and managers who seek to trace the sources of the formulae involved
in a computation.  A vignette is also provided for wider context.

A good way to start using `whalestrike` is to become familiar with `strike()`,
a function that simulates a collision event and produces an object that can
then be plotted (using the R generic function system) in multiple ways.

Three arguments must be provided to `strike()`.  The first sets the times at
which model output is desired, the second establishes the initial locations
and speeds of both ship and whale, and the third describes the biological and
physical properties of both ship and whale. The examples in the
documentation for `strike()` provide a good starting point for these three
arguments.  For example, the first example produced by typing
```R
library(whalestrike)
example(strike)
```
in an R console will run a sample simulation, and show a three-panel plot
(reproduced here as Figure 1) of three key aspects of the collision: the
relative position of whale and ship, the compression of the whale's layers, and
a lethality index developed by @kelley_assessing_2020.

![Diagram produced by typing `example(strike,package="whalestrike")` in an R session. Three of the possible 12 plots are shown. **Left:** positions of 45-tonne vessel initially moving at 10 knots (dashed line) and the boundaries between the three layers on the shipward side of the whale. Contact occurs at about $0.2$ s, and continues to about $0.6$ s Note that a small fishing vessel is used in this simulation, and so its speed is significantly reduced by the collision.  **Middle:** As the left panel, but showing only the positions of the interfaces between the whale layers, relative to the whale's centre position. In this view, the thin skin can be seen at the bottom, with blubber, sublayer and then bone to the interior. **Right:** An index of lethality, with the line thickened for the interval during which stresses are predicted to exceed a threshold for lethal damage according to @kelley_assessing_2020.](figure1.png)

Running the simulation and plotting the results as in Figure 1 takes a fraction
of a second on a three-year old laptop, showing that the system is
computationally inexpensive. This is helpful in detailed studies that involve
calling `strike()` with a wide suite of parameter values, such as the creation
of the diagrams in @kelley_assessing_2020, some of which involved tens of
thousands of model runs to cover parameter space in detail.

There also applications for which a few model runs may suffice. For such work,
`whalestrike` provides an R-shiny application that is run by executing
```R
library(whalestrike)
app()
```
in an R console. This creates a window (reproduced in Figure 2) in which there
are sliders and other GUI elements for controlling the simulation. This tool
makes it is easy to explore "what if" scenarios for ship strikes. For example,
adjusting the slider controlling ship speed and monitoring the lethality index
panel might be useful in motivating policies for a given ship class, and
altering the ship mass or prow geometry at a given speed casts some light on
the issue of whether one speed restriction ought to be applied to all types of
ship.


# Conclusions

The `whalestrike` package is intended to provide guidance for the development of
marine policies related to ship speeds.  It is written in R, a language that is
familiar to many marine biologists, and one that also offers a vast array of
statistical tools for analysing the results of model simulations.

Being founded on physical principles, `whalestrike` complements the more
statistical approaches of most studies in this field. An advantage of this
foundation is a reduction in the need for a database of collision events that
covers the relevant range of ship types and speeds.  This is especially
important for an application such as this, because each entry in a collision
database represents the loss of one more member of a species that is already at
high risk of extinction.

![View of an interactive application for simulaating ship-whale collisions. Normally some simple instructions appear at the top of this display, but these are trimmed for presentation here.  The labels on various GUI elements should give an indication of the properties that a user can vary with this tool.](figure2.png)


# Acknowledgements

I thank Sean W. Brillant for asking me to explain what I (then) knew about
ship-whale collisions, and James P. Vlasic for finding and collating published
records of such collisions and discussing representations of material
properties with me.  I would not have considered writing the `whalestrike`
package without the motivation of working with these two fine collaborators.
I am also grateful to Christopher T. Taggart, for discussions about whales and
other things that I never realized would interest me as much as they do.

# References
