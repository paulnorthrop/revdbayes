# Informative GEV prior on a quantile scale

Informative GEV prior for GEV parameters (\\\mu, \sigma, \xi\\)
constructed on the quantile scale. For information about how to set this
prior see
[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md).

## Usage

``` r
gev_quant(pars, prob, shape, scale, min_xi = -Inf, max_xi = Inf, trendsd = 0)
```

## Arguments

- pars:

  A numeric vector of length 3. GEV parameters (\\\mu, \sigma, \xi\\).

- prob:

  A numeric vector of length 3 containing exceedance probabilities
  (\\p_1, p_2, p_3\\) such that \\p_1 \> p_2 \> p_3\\. If the values in
  `quant` are not ordered from largest to smallest then they will be
  ordered inside `set_prior` without warning.

- shape, scale:

  Numeric vectors of length 3. Shape and scale parameters specifying
  (independent) gamma prior distributions placed on the differences
  between the quantiles corresponding to the probabilities given in
  `prob`.

- min_xi:

  A numeric scalar. Prior lower bound on \\\xi\\.

- max_xi:

  A numeric scalar. Prior upper bound on \\\xi\\.

- trendsd:

  Has no function other than to achieve compatibility with function in
  the evdbayes package.

## Value

The log of the prior density.

## Details

See Coles and Tawn (1996) and/or Stephenson (2016) for details.

Note that the lower end point of the distribution of the distribution of
the variable in question is assumed to be equal to zero. If this is not
the case then the user should shift the data to ensure that this is
true.

## References

Coles, S. G. and Tawn, J. A. (1996) A Bayesian analysis of extreme
rainfall data. *Appl. Statist.*, **45**, 463-478.

Stephenson, A. (2016) Bayesian inference for extreme value modelling. In
*Extreme Value Modeling and Risk Analysis: Methods and Applications*
(eds D. K. Dey and J. Yan), 257-280, Chapman and Hall, London.
[doi:10.1201/b19721](https://doi.org/10.1201/b19721) .
