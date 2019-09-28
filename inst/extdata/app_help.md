When this starts, the sliders and tick-boxes are set up to model a small
fishing boat, of mass 45 tonnes, moving at speed 10 knots towards a whale of
length 13.7m. (The whale length is used to compute its mass, using a formula
that is described by the output of typing
`help("whaleMassFromLength","whalestrike")` in an R console).

Sliders are provided for setting certain key properties of the ship and the
whale, with italic labels used for properties likely to be adjusted during
simulations.  The details of these and the other parameters are revealed by
typing `help("parameters","whalestrike")` and `help("strike","whalestrike")` in
an R console.

To the right of the sliders is a column of checkboxes that control the plotted
output. (For the details of the plots, type `help("plot.strike","plot")` in an
R console.) At startup, three of these boxes are ticked, yielding a display
with three panels showing the time history of the simulation.  The left-hand
plot panel shows whale and boat location, the former with an indication of the
interfaces between skin, blubber, sublayer, and bone. Peak ship and whale
accelerations are indicated with labels inside this panel.  The middle panel
shows the same information as the one to its left, but with a whale-centred
coordinate system, and with labels for the components. The right panel is an
indication of the estimated threat to the four layers of the whale, with curves
that are filled with colours that darken with the degree of threat. The
dividing lines are the quantiles of a logistic fit of published reports of
whale injury (with 0 meaning no injury or minor injury and 1 meaning severe or
fatal injury) to the base-10 logarithm of compressive stress.

Much can be learned by adjusting the sliders and examining the plotted output.
As an exercise, try setting to a particular ship mass of interest, and then to
slide the ship speed to higher and lower values, whilst monitoring the "threat"
panel. Next, try altering the blubber and sublayer thicknesses, e.g. using
smaller values to represent a strike at the whale mandible.  Having built some
intuition with these experiments, move on to altering the properties of the
ship, exploring the effect of changing ship speed, ship mass, and the width and
height of the impact zone.

