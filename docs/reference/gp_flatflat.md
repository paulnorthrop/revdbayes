# Flat prior for GP parameters (\\\sigma, \xi\\)

For information about this and other priors see
[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md).

## Usage

``` r
gp_flatflat(pars, min_xi = -Inf, max_xi = Inf, trendsd = 0, upper = NULL)
```

## Arguments

- pars:

  A numeric vector of length 2. GP parameters (\\\sigma, \xi\\).

- min_xi:

  A numeric scalar. Prior lower bound on \\\xi\\. Must not be `-Inf`
  because this results in an improper posterior.

- max_xi:

  A numeric scalar. Prior upper bound on \\\xi\\.

- trendsd:

  Has no function other than to achieve compatibility with function in
  the evdbayes package.

- upper:

  A positive numeric scalar. The upper endpoint of the GP distribution.

## Value

The log of the prior density.
