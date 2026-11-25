# Bivariate normal prior for GP parameters (\\log \sigma, \xi\\)

For information about this and other priors see
[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md).

## Usage

``` r
gp_norm(pars, mean, icov, min_xi = -Inf, max_xi = Inf, trendsd = 0)
```

## Arguments

- pars:

  A numeric vector of length 2. GP parameters (\\\sigma, \xi\\).

- mean:

  A numeric vector of length 2. Prior mean.

- icov:

  A 2x2 numeric matrix. The inverse of the prior covariance matrix.

- min_xi:

  A numeric scalar. Prior lower bound on \\\xi\\.

- max_xi:

  A numeric scalar. Prior upper bound on \\\xi\\.

- trendsd:

  Has no function other than to achieve compatibility with function in
  the evdbayes package.

## Value

The log of the prior density.
