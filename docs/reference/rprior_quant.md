# Prior simulation of GEV parameters - prior on quantile scale

Simulates from the prior distribution for GEV parameters proposed in
Coles and Tawn (1996), based on independent gamma priors for differences
between quantiles.

## Usage

``` r
rprior_quant(n, prob, shape, scale, lb = NULL, lb_prob = 0.001)
```

## Arguments

- n:

  A numeric scalar. The size of sample required.

- prob:

  A numeric vector of length 3. Exceedance probabilities corresponding
  to the quantiles used to specify the prior distribution. The values
  should *decrease* with the index of the vector. If not, the values in
  `prob` will be sorted into decreasing order without warning.

- shape:

  A numeric vector of length 3. Respective shape parameters of the gamma
  priors for the quantile differences.

- scale:

  A numeric vector of length 3. Respective scale parameters of the gamma
  priors for the quantile differences.

- lb:

  A numeric scalar. If this is not `NULL` then the simulation is
  constrained so that `lb` is an approximate lower bound on the GEV
  variable. Specifically, only simulated GEV parameter values for which
  the 100`lb_prob`% quantile is greater than `lb` are retained.

- lb_prob:

  A numeric scalar. The non-exceedance probability involved in the
  specification of `lb`. Must be in (0,1). If `lb=NULL` then `lb_prob`
  is not used.

## Value

An `n` by 3 numeric matrix.

## Details

The simulation is based on the way that the prior is constructed. See
Coles and Tawn (1996), Stephenson (2016) or the evdbayes user guide for
details of the construction of the prior. First, the quantile
differences are simulated from the specified gamma distributions. Then
the simulated quantiles are calculated. Then the GEV location, scale and
shape parameters that give these quantile values are found, by solving
numerically a set of three non-linear equations in which the GEV
quantile function evaluated at the values in `prob` is equated to the
simulated quantiles. This is reduced to a one-dimensional optimisation
over the GEV shape parameter.

## References

Coles, S. G. and Tawn, J. A. (1996) A Bayesian analysis of extreme
rainfall data. *Appl. Statist.*, **45**, 463-478.

Stephenson, A. 2016. Bayesian Inference for Extreme Value Modelling. In
*Extreme Value Modeling and Risk Analysis: Methods and Applications*,
edited by D. K. Dey and J. Yan, 257-80. London: Chapman and Hall.
[doi:10.1201/b19721](https://doi.org/10.1201/b19721)

## See also

[`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
and
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md)
for sampling from an extreme value posterior distribution.

## Examples

``` r
pr <- 10 ^ -(1:3)
sh <- c(38.9, 7.1, 47)
sc <- c(1.5, 6.3, 2.6)
x <- rprior_quant(n = 1000, prob = pr, shape = sh, scale = sc)
x <- rprior_quant(n = 1000, prob = pr, shape = sh, scale = sc, lb = 0)
```
