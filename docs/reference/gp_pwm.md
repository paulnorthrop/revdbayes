# Probability-weighted moments estimation of generalised Pareto parameters

Uses the methodology of Hosking and Wallis (1987) to estimate the
parameters of the generalised Pareto (GP) distribution.

## Usage

``` r
gp_pwm(gp_data, u = 0)
```

## Arguments

- gp_data:

  A numeric vector of raw data, assumed to be a random sample from a
  probability distribution.

- u:

  A numeric scalar. A threshold. The GP distribution is fitted to the
  excesses of `u`.

## Value

A list with components

- `est`: A numeric vector. PWM estimates of GP parameters \\\sigma\\
  (scale) and \\\xi\\ (shape).

- `se`: A numeric vector. Estimated standard errors of \\\sigma\\ and
  \\\xi\\.

- `cov`: A numeric matrix. Estimate covariance matrix of the the PWM
  estimators of \\\sigma\\ and \\\xi\\.

## References

Hosking, J. R. M. and Wallis, J. R. (1987) Parameter and Quantile
Estimation for the Generalized Pareto Distribution. Technometrics,
29(3), 339-349. [doi:10.2307/1269343](https://doi.org/10.2307/1269343) .

## See also

[`gp`](https://github.com/paulnorthrop/revdbayes/reference/gp.md) for
details of the parameterisation of the GP distribution.

## Examples

``` r
u <- quantile(gom, probs = 0.65)
gp_pwm(gom, u)
#> $est
#>      sigma         xi 
#> 1.77767265 0.08739445 
#> 
#> $se
#>     sigma        xi 
#> 0.2605594 0.1099028 
#> 
#> $cov
#>             sigma          xi
#> sigma  0.06789120 -0.02025464
#> xi    -0.02025464  0.01207863
#> 
```
