# GUI application for whale simulation

Graphical-user-interface tool for exploring whale-strike simulations.

## Usage

``` r
app(debug = FALSE)
```

## Arguments

- debug:

  logical value indicating whether to print output to the R console as
  the computation is done.

## Details

Sliders, buttons, and choosers are grouped into panes that appear on the
left of the view. When `app()` first opens, all of these panes are
closed. To get acquainted with the app, try adjusting the controllers
that *are* visible on the initial view. Then, open the "ship" pane and
increase the ship mass. Do you find that the results make qualitative
sense? Continue this process, exploring all the panes. A half-hour of
such exploration should be enough to build enough confidence to start
investigating practical applications. To learn more about how the
simulations are carried out, and to read more about the underlying goals
of this tool, please consult Kelley et al. (2021) and Kelley (2024).
Extensive details on the calculations are provided in the help pages for
the various functions of the whalestrike package, of which that for
[`whalestrike()`](https://dankelley.github.io/whalestrike/reference/whalestrike.md)
is a good starting point.

More information on `app()` in video form on
[youtube](https://youtu.be/kTMl3nXa5A4).

Note that an older version of a similar GUI application is still
available as
[`app_2025()`](https://dankelley.github.io/whalestrike/reference/app_2025.md),
but it is not maintained and is slated for removal in the early months
of 2026.

## References

Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. "Assessing the
Lethality of Ship Strikes on Whales Using Simple Biophysical Models."
Marine Mammal Science 37, no. 1 (2021): 251–67.
https://doi.org/10.1111/mms.12745.

Kelley, Dan E."“Whalestrike: An R Package for Simulating Ship Strikes on
Whales." Journal of Open Source Software 9, no. 97 (2024): 6473.
https://doi.org/10.21105/joss.06473.

Mayette, Alexandra. "Whale Layer Thickness." December 15, 2025.
(Personal communication of a 5-page document.)

Mayette, Alexandra, and Sean W. Brillant. "A Regression-Based Method to
Estimate Vessel Mass for Use in Whale-Ship Strike Risk Models." PloS One
21, no. 1 (2026): e0339760.
https://doi.org/10.1371/journal.pone.0339760.

## See also

Other interactive apps:
[`app_2025()`](https://dankelley.github.io/whalestrike/reference/app_2025.md)

## Author

Dan Kelley
