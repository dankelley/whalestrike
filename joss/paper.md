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
North Atlantic right whale. Although speed reduction policies have been
developed and implemented, simple considerations suggest that other factors,
such as ship mass and prow shape, should also be taken into account. An R
package called `whalestrike` has been developed to address such issues. It was
used by @kelley_assessing_2021 in the development of a biomechanically based
criterion for the lethality of ship strikes, but this was just a starting
point. The next step, and goal of the present paper, is to introduce the model
code to a broader community, encouraging its use and development by diverse
researchers and policy makers.

# Statement of need

Ship collisions pose significant threats to endangered marine mammals
[@laist_collisions_2001]. The case of the North Atlantic right whale
(*Eubalaena glacialis*) is of particular concern, since it is considered a
"Critically Endangered" species [@iucn_eubalaena_2020].  The trend is not
encouraging for this animal, the population of which was $340\pm7$ in 2021
[@pettis_north_2023], down significantly from $483$ in 2010
[@pace_state-space_2017].  According to necropsy studies, ship strikes are
involved in over half of right whale deaths [@campbell-malone_gross_2008].

Motivated by such studies, efforts have been made in recent years to mitigate
the consequences of collision by imposing speed restrictions on ships, and
evidence of successful results [e.g. @conn_vessel_2013] has led to marine
policy changes, including both static and dynamic zones of speed restriction
[e.g. @transport_canada_protecting_2023]. Even so, it seems unwise to measure
the success of speed restrictions by counting dead or injured whales, given the
low numbers alive today. In addition, basic reasoning reveals that speed cannot
be the only factor. Variations in ship mass, prow shape, etc. should also be
considered, along with whale morphometric parameters, making for a
multifactorial problem that requires even more data to achieve statistical
reliability.

With this in mind, @kelley_assessing_2021 devised a simple numerical model of
the biophysical dynamics of collisions between vessels and whales. This
involved using published records of whale injury and death (across multiple
species) to calibrate a lethality criterion that depends not just on vessel
speed, but also on other factors such as vessel mass, prow geometry, etc., in
addition to whale morphometric parameters. The model was expressed in the R
language, because it is familiar to many marine biologists and because it
provides a wide scope of statistical tools for followup work. The purpose of
the present paper is not to recapitulate the results of @kelley_assessing_2021,
but rather to introduce readers to the `whalestrike` code.

To the author's knowledge, there are no other R packages that address the topic
of collision mechanics and whale lethality in this way.  Readers interested in
the related topic of the probability of collision given shipping patterns and
whale distributions might find it helpful to start with R packages named
`whalemap` (@johnson_whalemap_2021) and `shipstrike`
(@keen_ericmkeenshipstrike_2023).

# Model formulation

A desire to produce a GUI (graphical user interface) tool permitting easy
exploration of various model scenarios led to a decision to create a simplified
model that yields results quickly, as opposed to a much more
computationally-expensive finite element model that might account more
accurately for the deformation of whale flesh [see e.g.
@raymond_development_2007].  The new model ignores ship deformation upon
impact, and considers whale deformation to occur only in a specified impact
area dictated by the geometry of the ship's prow.  It uses a layered scheme to
represent the whale, with skin covering blubber, with that blubber covering a
sub-layer representing a combination of muscle and organs, and with bone at the
core. Thicknesses and material properties (including a nonlinear stress-strain
relationship) for each layer are taken from the literature, with the idea being
that adjusting these parameters will provide a way to simulate strikes on
different species, or at different body locations. Skin deformation is modeled
with both extension and compression forces, while only compression is
considered for interior layers. Lacking reliable information on failure limits
for these biomaterials, critical values for stresses are posited in the context
    of published results of the damage experienced in documented ship strikes.
    @kelley_assessing_2021 provide more information on the methodologies
    and issues involved.

The model is simple, with acceleration and forces being linked via Newton's
second law. Numerical integration is carried out with the `lsoda()` function of
the `deSolve` R package [@soetaert_solving_2010]. Predicted velocities are used
in computing water drag and the relative positions of ship and whale. After
contact is made, extensive and compressive forces arise, and the stresses
associated with these forces are monitored in the context of the lethality
index proposed by @kelley_assessing_2021. Dynamical elements involved in this
process are distilled into several dozen model parameters. Although a great deal of
effort has been put into formulating these parameters appropriately (mainly for
application to right whales), `whalestrike` offers a simple way for users to
adjust each of them, if needed for new applications.

# Package installation and use

The package may be installed from its development website by typing
```R
library(remotes)
install_github("dankelley/whalestrike", ref="main")
```
in an R console^[There are plans to submit `whalestrike` to the Comprehensive R
Archive Network [@cran_comprehensive_2024].].  Once this is done, the user has
access to over 20 functions and their documentation.  The latter makes frequent
reference to the scientific literature, since the package is aimed at
scientists and managers who seek to trace the sources of the formulae involved
in a computation.  A vignette is also provided for wider context.

A good way to start using `whalestrike` is to become familiar with `strike()`,
a function that simulates a collision event and produces an object that can
then be plotted (using the R generic function system) in multiple ways.

Three arguments must be provided to `strike()`.  The first sets the times at
which model output is desired, the second establishes the initial locations and
speeds of the ship and the whale, and the third describes the biological and
physical properties of the ship and the whale. An example is provided
in the documentation of `strike()`, e.g. typing


```R
library(whalestrike)
example(strike)
```

in an R console will run a sample simulation, and show a three-panel plot
(reproduced here as Figure 1) of the relative position of whale and ship, the
compression of the whale's layers, and the associated lethality index.

![Diagram produced by typing `example(strike,package="whalestrike")` in an R session. Three of the possible 12 plots are shown. **Left:** positions of 45-tonne vessel initially moving at 10 knots (dashed line) and the boundaries between the three layers on the shipward side of the whale. Contact occurs at about $0.2$ s, and continues to about $0.6$ s. Note that a small fishing vessel is used in this simulation, and so its speed is significantly reduced by the collision.  **Middle:** As the left panel, but showing only the positions of the interfaces between the whale layers, relative to the whale's centre position. In this view, the thin skin can be seen at the bottom, with blubber, sublayer and then bone to the interior. **Right:** An index of lethality, with the line thickened for the interval during which stresses are predicted to exceed a threshold for lethal damage according to @kelley_assessing_2021.](figure1.png)

Running the simulation and plotting the results as in Figure 1 takes a fraction
of a second on a typical laptop, showing that the system is computationally
inexpensive. This is helpful in detailed studies that involve calling
`strike()` with a wide suite of parameter values, such as the creation of the
diagrams in @kelley_assessing_2021, some of which involved of order $10^5$
model runs to cover parameter space in detail.

There are also applications for which a few model runs may suffice. For such
work, `whalestrike` provides an R-shiny application that is run by executing
the following in an R console.

```R
library(whalestrike)
app2()
```

Doing this creates a window, shown in Figure 2 here and in video form online
[@dan_kelley_using_2024]. With a variety of sliders and other GUI elements, the
application makes it easy to explore "what if" scenarios for ship strikes. For
example, monitoring the Lethality Index plot while adjusting the
ship-specification tools and the impact speed might prove useful in discussions
of speed restrictions across ship masses and classes. More detailed (and
reproducible) work is also facilitated by `app2()`, because it can display the
underlying code used in the simulation, providing a good starting point for
more extensive analyses, such as the exploration of covarying parameters.

![View of an interactive application for simulating ship-whale collisions. Adjusting the slider for ship speed will reveal that the vessel in this simulation would have to slow down to 6.5 knots in order to reduce the inferred Lethality Index below a critical value (dashed line in right panel).](figure2.png)


# Conclusions

The `whalestrike` package is intended to provide guidance for the development
of marine policies related to ship speeds.  It is written in R, a language that
is familiar to many marine biologists, and one that also offers a vast array of
statistical tools for analysing the results of model simulations. In addition
to tools for detailed control of simulations, the package also provides a
GUI-based tools for basic exploration, which may be useful in making or
explaining policy decisions.

Being founded on physical principles, `whalestrike` complements the more
statistical approaches of most studies in this field. An advantage of this
foundation is a reduction in the need for a database of collision events that
covers the relevant range of ship types and speeds.  This is especially
important for an application such as this, because each entry in a collision
database represents the loss of one more member of a species that is already at
high risk of extinction.

# Acknowledgements

I thank Sean W. Brillant for asking me to explain what I (then) knew about
ship-whale collisions, and James P. Vlasic for finding and collating published
records of such collisions and discussing representations of material
properties with me. I would not have considered writing the `whalestrike`
package without the motivation of working with these two fine collaborators.
Jaimie Harbin provided useful comments on this manuscript.  I am also grateful
to Christopher T. Taggart, for discussions about whales and other things that I
never realized would interest me as much as they do. Finally, I thank the
reviewers of this manuscript, for careful reading of this document and
thoughtful examination of the `whalestrike` package code.

# References
