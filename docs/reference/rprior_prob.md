# Prior simulation of GEV parameters - prior on probability scale

Simulates from the prior distribution for GEV parameters based on
Crowder (1992), in which independent beta priors are specified for
ratios of probabilities (which is equivalent to a Dirichlet prior on
differences between these probabilities).

## Usage

``` r
rprior_prob(n, quant, alpha, exc = FALSE, lb = NULL, lb_prob = 0.001)
```

## Arguments

- n:

  A numeric scalar. The size of sample required.

- quant:

  A numeric vector of length 3. Contains quantiles \\q_1, q_2, q_3\\. A
  prior distribution is placed on the non-exceedance (`exc = FALSE`) or
  exceedance (`exc = TRUE`) probabilities corresponding to these
  quantiles. The values should *increase* with the index of the vector.
  If not, the values in `quant` will be sorted into increasing order
  without warning.

- alpha:

  A numeric vector of length 4. Parameters of the Dirichlet distribution
  for the exceedance probabilities.

- exc:

  A logical scalar. Let \\M\\ be the GEV variable, \\r_q = P(M \leq
  q)\\, \\p_q = P(M \> q) = 1 - r_q\\ and `quant` = (\\q_1, q_2, q_3\\).
  If `exc = FALSE` then a Dirichlet(`alpha`) distribution is placed on
  \\(r\_{q_1}, r\_{q_2} - r\_{q_1}, r\_{q_3} - r\_{q_2}, 1 -
  r\_{q_3})\\, as in Northrop et al. (2017). If `exc = TRUE` then a
  Dirichlet(`alpha`) distribution is placed on \\(1 - p\_{q_1},
  p\_{q_1} - p\_{q_2}, p\_{q_2} - p\_{q_3}, p\_{q_3})\\, where \\p_q =
  P(M \> q)\\, as in Stephenson (2016).

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
Stephenson (1996) the evdbayes user guide or Northrop et al. (2017)
Northrop et al. (2017) for details of the construction of the prior.
First, differences between probabilities are simulated from a Dirichlet
distribution. Then the GEV location, scale and shape parameters that
correspond to these quantile values are found, by solving numerically a
set of three non-linear equations in which the GEV quantile function
evaluated at the simulated probabilities is equated to the quantiles in
`quant`. This is reduced to a one-dimensional optimisation over the GEV
shape parameter.

## References

Crowder, M. (1992) Bayesian priors based on parameter transformation
using the distribution function. *Ann. Inst. Statist. Math.*, **44**(3),
405-416. <https://link.springer.com/article/10.1007/BF00050695>

Stephenson, A. 2016. Bayesian Inference for Extreme Value Modelling. In
*Extreme Value Modeling and Risk Analysis: Methods and Applications*,
edited by D. K. Dey and J. Yan, 257-80. London: Chapman and Hall.
[doi:10.1201/b19721](https://doi.org/10.1201/b19721)

Northrop, P. J., Attalides, N. and Jonathan, P. (2017) Cross-validatory
extreme value threshold selection and uncertainty with application to
ocean storm severity. *Journal of the Royal Statistical Society Series
C: Applied Statistics*, **66**(1), 93-120.
[doi:10.1111/rssc.12159](https://doi.org/10.1111/rssc.12159)

## See also

[`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
and
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md)
for sampling from an extreme value posterior distribution.

## Examples

``` r
quant <- c(85, 88, 95)
alpha <- c(4, 2.5, 2.25, 0.25)
x <- rprior_prob(n = 1000, quant = quant, alpha = alpha, exc = TRUE)
x <- rprior_prob(n = 1000, quant = quant, alpha = alpha, exc = TRUE, lb = 0)
```
