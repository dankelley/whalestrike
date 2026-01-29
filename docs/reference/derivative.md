# Calculate derivative using first difference

The derivative is estimated as the ratio of the first-difference of
`var` divided by the first-difference of `time`. To make the results
have the same length as `time`, the final result is appended at the end.

## Usage

``` r
derivative(var, t)
```

## Arguments

- var:

  variable.

- t:

  time in seconds.

## Value

Derivative estimated by using
[`diff()`](https://rdrr.io/r/base/diff.html) on both `var` and `time`.

## Author

Dan Kelley
