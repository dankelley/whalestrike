# Summarize a parameters object

This provides an overview of the contents of an object created with
[`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md).

## Usage

``` r
# S3 method for class 'parameters'
summary(object, ...)
```

## Arguments

- object:

  an object of class `"parameters"`, as created with
  [`parameters()`](https://dankelley.github.io/whalestrike/reference/parameters.md).

- ...:

  ignored

## Author

Dan Kelley

## Examples

``` r
summary(parameters())
#> Whale and ship properties, as created by parameters():
#>   Ship properties:
#>     ms:          45000 kg  -- mass
#>     Ss:         154.54 m^2 -- wetted area
#>     Ly:           1.15 m   -- width of impact area
#>     Lz:           1.15 m   -- height of impact area
#>     Cs:           0.01     -- drag coefficient
#>   Whale properties:
#>     mw:        29993.9 kg  -- mass
#>     Sw:        64.6723 m^2 -- wetted area
#>     l[1]:        0.025 m   -- thickness of skin
#>     l[2]:         0.16 m   -- thickness of blubber
#>     l[3]:         1.12 m   -- thickness of sublayer
#>     l[4]:          0.1 m   -- half-thickness of bone
#>     lsum:        1.405 m   -- a[1]+a[2]+a[3]+a[4]
#>     a[1]:    1.780e+08 Pa  -- skin stress factor
#>     a[2]:    1.580e+05 Pa  -- blubber stress factor
#>     a[3]:    1.580e+05 Pa  -- sublayer stress factor
#>     a[4]:    8.540e+09 Pa  -- bone stress factor
#>     b[1]:          0.1     -- skin stress nonlinearity term
#>     b[2]:         2.54     -- blubber stress nonlinearity term
#>     b[3]:         2.54     -- sublayer stress nonlinearity term
#>     b[4]:          0.1     -- bone stress nonlinearity term
#>     s[1]:    1.960e+07 Pa  -- skin strength
#>     s[2]:    2.550e+05 Pa  -- blubber strength
#>     s[3]:    2.550e+05 Pa  -- sublayer strength
#>     s[4]:    2.290e+07 Pa  -- bone strength
#>     theta:          55 deg -- impact dimple angle
#>     Cw:       0.0025       -- drag coefficient
#>   Functions:
#>     stressFromStrain()     -- function to compute stress
#>     logistic()             -- function to compute lethality index
```
