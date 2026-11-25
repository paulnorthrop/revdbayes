# Maximal data information (MDI) prior for GEV parameters (\\\mu, \sigma, \xi\\)

For information about this and other priors see
[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md).

## Usage

``` r
gev_mdi(pars, a = 0.577215664901532, min_xi = -1, max_xi = Inf, trendsd = 0)
```

## Arguments

- pars:

  A numeric vector of length 3. GEV parameters (\\\mu, \sigma, \xi\\).

- a:

  A numeric scalar. The default value, Euler's constant, gives the MDI
  prior.

- min_xi:

  A numeric scalar. Prior lower bound on \\\xi\\. Must not be `-Inf`
  because this results in an improper posterior.

- max_xi:

  A numeric scalar. Prior upper bound on \\\xi\\.

- trendsd:

  Has no function other than to achieve compatibility with function in
  the evdbayes package.

## Value

The log of the prior density.
