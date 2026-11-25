# Maximum likelihood estimation of generalised Pareto parameters

Uses the methodology of Grimshaw (1993) to find the MLEs of the
parameters of the generalised Pareto distribution, based on a sample of
positive values. The function is essentially the same as that made
available with Grimshaw (1993), with only minor modifications.

## Usage

``` r
grimshaw_gp_mle(x)
```

## Arguments

- x:

  A numeric vector containing only **positive** values, assumed to be a
  random sample from a generalized Pareto distribution.

## Value

A numeric vector of length 2. The estimates of the **negated** shape
parameter \\k (= -\xi)\\ and the scale parameter \\a (= \sigma)\\.

## References

Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates for the
Generalized Pareto Distribution. Technometrics, 35(2), 185-191. and
Computing (1991) 1, 129-133.
[doi:10.1080/00401706.1993.10485040](https://doi.org/10.1080/00401706.1993.10485040)
.

## See also

[`gp`](https://github.com/paulnorthrop/revdbayes/reference/gp.md) for
details of the parameterisation of the GP distribution, in terms of
\\\sigma\\ and \\\xi\\.

## Examples

``` r
u <- quantile(gom, probs = 0.65)
grimshaw_gp_mle((gom - u)[gom > u])
#> $k
#> [1] -0.101175
#> 
#> $a
#> [1] 1.751086
#> 
```
