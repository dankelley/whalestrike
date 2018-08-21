When this starts, the sliders and tick-boxes are set up to model a small
fishing boat, of mass 45 tonnes, moving at speed 10 knots towards a whale of
length 13.7m (used with a corrected version of the formula for North Atlantic
Right Whales provided by Fortune et al, 2012).

Sliders are provided for setting certain key properties of the ship and the
whale, with red labels used for those properties that are deemed most likely to
be adjusted during simulations.  The details of these and the other parameters
can be found by typing `help("parameters","whalestrike")` and
`help("strike","whalestrike")` in the R console.

To the right of the sliders is a column of checkboxes that control the plotted
output. At startup, three of these boxes are ticked, yielding a display with
three panels showing the time history of the simulation.  The left-hand plot
panel shows whale and boat location, the former with an indication of the
interfaces between skin, blubber, sublayer, and bone. The middle panel shows
the same information as the left one, but with a whale-centred coordinate
system, and with labels for the components. The right panel is an indication of
the estimated threat to the four layers of the whale, with curves that are
filled with grey for time intervals when the impact stress (force/area) is less
than the strength of the material in the layer, and black for times when that
threshold is exceeded. The details of these and the other possible panels can
be found by typing `help("plot.strike","whalestrike")` in the R console.

Much can be learned by adjusting the coloured sliders and examining the plotted
output. As an exercise, try setting to a particular ship mass of interest, and
then to slide the ship speed to higher and lower values, whilst monitoring the
"threat" panel for black regions. This will reveal a critical speed for
conditions that threaten the whale.  Next, try altering the sublayer thickness,
which is a surrogate for location along the whale body, because e.g. the
sublayer is thinner near the mandible.

Advanced users are likely to want to alter the values of impact width and
height. The default setting are intended to mimic a small fishing boat, such as
a Cape Islander. Try lowering the width, to simulate a strike by a daggerboard
or keel of a sailing boat.

**References**

* Fortune, Sarah M. E., Andrew W. Trites, Wayne L. Perryman, Michael J. Moore,
Heather M. Pettis, and Morgan S. Lynn. “Growth and Rapid Early Development of
North Atlantic Right Whales (Eubalaena Glacialis).” Journal of Mammalogy 93,
no. 5 (2012): 1342–54. https://doi.org/10.1644/11-MAMM-A-297.1.

