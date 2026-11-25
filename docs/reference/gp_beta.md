# Beta-type prior for GP shape parameter \\\xi\\

For information about this and other priors see
[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md).

## Usage

``` r
gp_beta(pars, min_xi = -1/2, max_xi = 1/2, pq = c(6, 9), trendsd = 0)
```

## Arguments

- pars:

  A numeric vector of length 2. GP parameters (\\\sigma, \xi\\).

- min_xi:

  A numeric scalar. Prior lower bound on \\\xi\\.

- max_xi:

  A numeric scalar. Prior upper bound on \\\xi\\.

- pq:

  A numeric vector of length 2. See
  [`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md)
  for details.

- trendsd:

  Has no function other than to achieve compatibility with function in
  the evdbayes package.

## Value

The log of the prior density.
