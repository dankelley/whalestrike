# Pin numerical values between stated limits

Pin numerical values between stated limits

## Usage

``` r
pin(x, lower = NULL, upper = NULL)
```

## Arguments

- x:

  Vector or matrix of numerical values

- lower:

  Numerical values of minimum value allowed; set to `NULL` to avoid
  trimming the lower limit.

- upper:

  As for `lower`, but for the upper limit.

## Value

Copy of `x`, with any value that exceeds `lim` having been replaced by
`lim`.

## Author

Dan Kelley
