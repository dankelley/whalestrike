# whalestrike: A Package to Simulate Ship-Whale Collisions

This package solves Newton's second law for a simple model of a ship
colliding with a whale. This is a stripped-down model that does not
attempt to simulate the biomechanical interactions that can be simulated
in finite-element treatments such as that of Raymond (2007). For an
in-depth discussion of the reason for writing the model, of the
principles involved in its framing, and its use in developing a
criterion for strike lethality, please see Kelley et al. (2021).

## Details

The goal of the model is to establish a convenient framework for rapid
computation of impacts in a wide variety of conditions. The model runs
quickly enough to keep up with mouse movements to select conditions, in
a companion R shiny app. That app lets the user see the effects of
changing contact area, ship speed, etc., as a way to build intuition for
scenarios ranging from collision with a slow-moving fishing boat to
collision with a thin dagger-board of a much swifter racing sailboat.
Another advantage of the simple formulation is that it makes it easy to
modify various dynamical and biomechanical parameters, to add new
forces, and to explore a range of criteria for whale damage.

The documentation for
[`strike()`](https://dankelley.github.io/whalestrike/reference/strike.md)
provides a practical example of using the main functions of this
package, while the package vignette provides a general overview. Kelley
et al (2021) provide more detail about the mathematical framework of the
package, along with a discussion of its purpose and application to
real-world problems of ship strikes on whales.

## Further reading

- Daoust, Pierre-Yves, Emilie L. Couture, Tonya Wimmer, and Laura
  Bourque. "Incident Report. North Atlantic Right Whale Mortality Event
  in the Gulf of St. Lawrence, 2017." Canadian Wildlife Health
  Cooperative, Marine Animal Response Society, and Fisheries and Oceans
  Canada, 2018.
  <https://publications.gc.ca/site/eng/9.850838/publication.html>.

- Fortune, Sarah M. E., Andrew W. Trites, Wayne L. Perryman, Michael J.
  Moore, Heather M. Pettis, and Morgan S. Lynn. "Growth and Rapid Early
  Development of North Atlantic Right Whales (Eubalaena Glacialis)."
  Journal of Mammalogy 93, no. 5 (2012): 1342-54.
  [doi:10.1644/11-MAMM-A-297.1](https://doi.org/10.1644/11-MAMM-A-297.1)
  .

- Grear, Molly E., Michael R. Motley, Stephanie B. Crofts, Amanda E.
  Witt, Adam P. Summers, and Petra Ditsche. "Mechanical Properties of
  Harbor Seal Skin and Blubber–a Test of Anisotropy." Zoology 126
  (2018): 137-44.
  [doi:10.1016/j.zool.2017.11.002](https://doi.org/10.1016/j.zool.2017.11.002)
  .

- Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. "Assessing the
  Lethality of Ship Strikes on Whales Using Simple Biophysical Models."
  Marine Mammal Science 37, no. 1 (January 2021): 251–67.
  [doi:10.1111/mms.12745](https://doi.org/10.1111/mms.12745) .

- Kelley, Dan E. "Composite Spring," May 28, 2018.
  20180528_composite_string. Dan Kelley's working notes.

- Kelley, Dan. "Whale Area," June 23, 2018. 20180623_whale_area. Dan
  Kelley's working notes.

- Kelley, Dan. "Ship Propulsion," July 1, 2018.
  20180701_ship_propulsion. Dan Kelley's working notes.

- Kelley, Dan. "Whale Mass," July 7, 2018. 20180707_whale_mass. Dan
  Kelley's working notes.

- Kelley, Dan E."“Whalestrike: An R Package for Simulating Ship Strikes
  on Whales." Journal of Open Source Software 9, no. 97 (2024): 6473.
  https://doi.org/10.21105/joss.06473.

- MAN Diesel & Turbo. "Basic Principles of Propulsion." MAN Diesel &
  Turbo, 2011.
  `https://spain.mandieselturbo.com/docs/librariesprovider10/sistemas-propulsivos-marinos/basic-principles-of-ship-propulsion.pdf?sfvrsn=2`

- Manen, J. D. van, and P. van Oossanen. "Resistance." In Principles of
  Naval Architecture (Second Revision), Volume II - Resistance,
  Propulsion and Vibration, edited by Edward V Lewis, Second Edition,
  1-125. Jersey City, NJ: Society of Naval Architects and Marine
  Engineers (U.S.), 1988.

- Mayette, Alexandra. "Whale Layer Thickness." December 15, 2025.
  (Personal communication of a 5-page document.)

- Mayette, Alexandra, and Sean W. Brillant. "A Regression-Based Method
  to Estimate Vessel Mass for Use in Whale-Ship Strike Risk Models."
  PloS One 21, no. 1 (2026): e0339760.
  https://doi.org/10.1371/journal.pone.0339760.

- Miller, Carolyn A., Desray Reeb, Peter B. Best, Amy R. Knowlton,
  Moira W. Brown, and Michael J. Moore. "Blubber Thickness in Right
  Whales Eubalaena Glacialis and Eubalaena Australis Related with
  Reproduction, Life History Status and Prey Abundance." Marine Ecology
  Progress Series 438 (2011): 267-83.

- Moore, M.J., A.R. Knowlton, S.D. Kraus, W.A. McLellan, and R.K. Bonde.
  "Morphometry, Gross Morphology and Available Histopathology in North
  Atlantic Right Whale (Eubalaena Glacialis) Mortalities (1970 to
  2002)." Journal of Cetacean Research and Management 6, no. 3 (2005):
  199-214.

- Ng, Laurel J., Vladislav Volman, Melissa M. Gibbons, Pi Phohomsiri,
  Jianxia Cui, Darrell J. Swenson, and James H. Stuhmiller. "A
  Mechanistic End-to-End Concussion Model That Translates Head
  Kinematics to Neurologic Injury." Frontiers in Neurology 8, no. JUN
  (2017): 1-18.
  [doi:10.3389/fneur.2017.00269](https://doi.org/10.3389/fneur.2017.00269)

- Raymond, J. J. "Development of a Numerical Model to Predict Impact
  Forces on a North Atlantic Right Whale during Collision with a
  Vessel." University of New Hampshire, 2007.
  <https://scholars.unh.edu/thesis/309/>.

- Soetaert, Karline, Thomas Petzoldt, and R. Woodrow Setzer. "Solving
  Differential Equations in R: Package DeSolve." Journal of Statistical
  Software; Vol 1, Issue 9, 2010.
  [doi:10.18637/jss.v033.i09](https://doi.org/10.18637/jss.v033.i09) .

## See also

Useful links:

- <https://dankelley.github.io/whalestrike/>

- Report bugs at <https://github.com/dankelley/whalestrike/issues>

## Author

**Maintainer**: Dan Kelley <dan.kelley@dal.ca>
([ORCID](https://orcid.org/0000-0001-7808-5911))

Other contributors:

- James Vlasic <jvlasic@dal.ca>
  ([ORCID](https://orcid.org/0000-0002-3846-4391)) \[research team
  member\]

- Sean Brilliant <seanb@cwf-fcf.org>
  ([ORCID](https://orcid.org/0000-0001-5494-3475)) \[research team
  member\]

- Alexandra Mayette <alexandram@cwf-fcf.org>
  ([ORCID](https://orcid.org/0000-0002-9766-9565)) \[research team
  member\]
