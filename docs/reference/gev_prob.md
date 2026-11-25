# Informative GEV prior on a probability scale

Constructs an informative prior for GEV parameters (\\\mu, \sigma,
\xi\\), constructed on the probability scale. For information about how
to set this prior see
[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md).

## Usage

``` r
gev_prob(pars, quant, alpha, min_xi = -Inf, max_xi = Inf, trendsd = 0)
```

## Arguments

- pars:

  A numeric vector of length 3. GEV parameters (\\\mu, \sigma, \xi\\).

- quant:

  A numeric vector of length 3 containing quantiles (\\q_1, q_2, q_3\\)
  such that \\q_1 \< q_2 \< q_3\\. If the values in `quant` are not
  ordered from smallest to largest then they will be ordered inside
  `set_prior` without warning.

- alpha:

  A numeric vector of length 4. Parameters specifying a prior
  distribution for probabilities related to the quantiles in `quant`.
  See **Details** below.

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

A prior for GEV parameters \\(\mu, \sigma, \xi)\\, based on Crowder
(1992). This construction is typically used to set an informative prior,
based on specified quantiles \\q_1, q_2, q_3\\. There are two
interpretations of the parameter vector `alpha` = \\(\alpha_1, \alpha_2,
\alpha_3, \alpha_4)\\: as the parameters of beta distributions for ratio
of exceedance probabilities (Stephenson, 2016) and as the parameters of
a Dirichlet distribution for differences between non-exceedance
probabilities (Northrop et al., 2017). See these publications for
details.

## References

Crowder, M. (1992) Bayesian priors based on parameter transformation
using the distribution function *Ann. Inst. Statist. Math.*, **44**,
405-416. <https://link.springer.com/article/10.1007/BF00050695>.

Northrop, P. J., Attalides, N. and Jonathan, P. (2017) Cross-validatory
extreme value threshold selection and uncertainty with application to
ocean storm severity. *Journal of the Royal Statistical Society Series
C: Applied Statistics*, **66**(1), 93-120.
[doi:10.1111/rssc.12159](https://doi.org/10.1111/rssc.12159)

Stephenson, A. (2016) Bayesian inference for extreme value modelling. In
*Extreme Value Modeling and Risk Analysis: Methods and Applications*
(eds D. K. Dey and J. Yan), 257-280, Chapman and Hall, London.
[doi:10.1201/b19721](https://doi.org/10.1201/b19721) .

## See also

[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md)
for setting a prior distribution.

[`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
and
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md)
for sampling from an extreme value posterior distribution.

Sets the same prior as the function `prior.prob` in the evdbayes
package.
