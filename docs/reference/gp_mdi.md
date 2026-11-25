# Maximal data information (MDI) prior for GP parameters (\\\sigma, \xi\\)

For information about this and other priors see
[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md).

## Usage

``` r
gp_mdi(pars, a = 1, min_xi = -1, max_xi = Inf, trendsd = 0)
```

## Arguments

- pars:

  A numeric vector of length 2. GP parameters (\\\sigma, \xi\\).

- a:

  A numeric scalar. The default value of 1 gives the MDI prior.

- min_xi:

  A numeric scalar. Prior lower bound on \\\xi\\. Must not be `-Inf`
  because this results in an improper posterior. See Northrop and
  Attalides (2016) for details.

- max_xi:

  A numeric scalar. Prior upper bound on \\\xi\\.

- trendsd:

  Has no function other than to achieve compatibility with function in
  the evdbayes package.

## Value

The log of the prior density.

## References

Northrop, P.J. and Attalides, N. (2016) Posterior propriety in Bayesian
extreme value analyses using reference priors *Statistica Sinica*,
**26(2)**, 721–743
[doi:10.5705/ss.2014.034](https://doi.org/10.5705/ss.2014.034) .
