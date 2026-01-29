# GUI application for whale simulation

This is a replacement to the older
[`app_2025()`](https://dankelley.github.io/whalestrike/reference/app_2025.md),
which had a more awkward interface, and which is no longer maintained.

## Usage

``` r
app(debug = FALSE)
```

## Arguments

- debug:

  logical value indicating whether to print output to the R console as
  the computation is done.

## Details

Compared with
[`app_2025()`](https://dankelley.github.io/whalestrike/reference/app_2025.md),
the present function lacks the ability to save settings and reload them
later. This is mainly because it only works with locally-run operations,
not from server-run operations. The latter would require extra coding to
set up user's storage space, to prevent against web attacks, etc., which
is beyond the present purpose. However, there is an addition with
`app()` that might prove more useful: a button to display the code
required to reproduce the simulated state. This may be of help to the
those seeking to explore the results of simulations more precisely and
with greater reproducibility.

Sliders, buttons, and choosers are grouped into panes that appear on the
left of the view. When `app()` first opens, all of these panes are
closed. To get acquainted with the app, try adjusting the controllers
that *are* visible on the initial view. Then, open the "ship" pane and
increase the ship mass. Do you find that the results make qualitative
sense? Continue this process, exploring all the panes. It is hoped that
a half hour of such exploration will let users start to investigate
practical applications. For more about how the simulations are carried
out, as well as comments on some applications that may be of interest,
please consult Kelley et al. (2021).

More information on `app()` in video form on
[youtube](https://youtu.be/kTMl3nXa5A4).

## References

Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. "Assessing the
Lethality of Ship Strikes on Whales Using Simple Biophysical Models."
Marine Mammal Science, 37(1), 2021.
[doi:10.1111/mms.12745](https://doi.org/10.1111/mms.12745) .

## See also

Other interactive apps:
[`app_2025()`](https://dankelley.github.io/whalestrike/reference/app_2025.md)

## Author

Dan Kelley
