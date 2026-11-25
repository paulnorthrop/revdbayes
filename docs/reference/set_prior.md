# Construction of prior distributions for extreme value model parameters

Constructs a prior distribution for use as the argument `prior` in
[`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
and
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md).
The user can either specify their own prior function, returning the
**log** of the prior density, (using an R function or an external
pointer to a compiled C++ function) and arguments for hyperparameters or
choose from a list of in-built model-specific priors. Note that the
arguments `model = "gev"`, `model = "pp"` and `model = "os"` are
equivalent because a prior is specified is the GEV parameterisation in
each of these cases. Note also that for `model = "pp"` the prior GEV
parameterisation relates to the value of `noy` subsequently supplied to
[`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
or
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md).
The argument `model` is used for consistency with `rpost`.

## Usage

``` r
set_prior(
  prior = c("norm", "loglognorm", "mdi", "flat", "flatflat", "jeffreys", "beta", "prob",
    "quant"),
  model = c("gev", "gp", "pp", "os"),
  ...
)
```

## Arguments

- prior:

  Either

  - An R function, or a pointer to a user-supplied compiled C++
    function, that returns the value of the log of the prior density
    (see **Examples**), or

  - A character string giving the name of the prior. See **Details** for
    a list of priors available for each model.

- model:

  A character string. If `prior` is a character string then `model`
  gives the extreme value model to be used. Using either
  `model = "gev"`, `model = "pp"` or `model = "os"` will result in the
  same (GEV) parameterisation. If `prior` is a function then the value
  of `model` is stored so that in the subsequent call to `rpost`,
  consistency of the prior and extreme value model parameterisations can
  be checked.

- ...:

  Further arguments to be passed to the user-supplied or in-built prior
  function. For details of the latter see **Details** and/or the
  relevant underlying function:
  [`gp_norm`](https://github.com/paulnorthrop/revdbayes/reference/gp_norm.md),
  [`gp_mdi`](https://github.com/paulnorthrop/revdbayes/reference/gp_mdi.md),
  [`gp_flat`](https://github.com/paulnorthrop/revdbayes/reference/gp_flat.md),
  [`gp_flatflat`](https://github.com/paulnorthrop/revdbayes/reference/gp_flatflat.md),
  [`gp_jeffreys`](https://github.com/paulnorthrop/revdbayes/reference/gp_jeffreys.md),
  [`gp_beta`](https://github.com/paulnorthrop/revdbayes/reference/gp_beta.md),
  [`gev_norm`](https://github.com/paulnorthrop/revdbayes/reference/gev_norm.md),
  [`gev_loglognorm`](https://github.com/paulnorthrop/revdbayes/reference/gev_loglognorm.md),
  [`gev_mdi`](https://github.com/paulnorthrop/revdbayes/reference/gev_mdi.md),
  [`gev_flat`](https://github.com/paulnorthrop/revdbayes/reference/gev_flat.md),
  [`gev_flatflat`](https://github.com/paulnorthrop/revdbayes/reference/gev_flatflat.md),
  [`gev_beta`](https://github.com/paulnorthrop/revdbayes/reference/gev_beta.md),
  [`gev_prob`](https://github.com/paulnorthrop/revdbayes/reference/gev_prob.md),
  [`gev_quant`](https://github.com/paulnorthrop/revdbayes/reference/gev_quant.md).
  All these priors have the arguments `min_xi` (prior lower bound on
  \\\xi\\) and `max_xi` (prior upper bound on \\\xi\\).

## Value

A list with class `"evprior"`. The first component is the input prior,
i.e. either the name of the prior or a user-supplied function. The
remaining components contain the numeric values of any hyperparameters
in the prior.

## Details

Of the in-built named priors available in revdbayes only those specified
using `prior = "prob"`
([`gev_prob`](https://github.com/paulnorthrop/revdbayes/reference/gev_prob.md)),
`prior = "quant"`
([`gev_quant`](https://github.com/paulnorthrop/revdbayes/reference/gev_quant.md))
`prior = "norm"`
([`gev_norm`](https://github.com/paulnorthrop/revdbayes/reference/gev_norm.md))
or `prior = "loglognorm"`
([`gev_loglognorm`](https://github.com/paulnorthrop/revdbayes/reference/gev_loglognorm.md))
are proper. If `model = "gev"` these priors are equivalent to priors
available in the evdbayes package, namely `prior.prob`, `prior.quant`,
`prior.norm` and `prior.loglognorm`.

The other in-built priors are improper, that is, the integral of the
prior function over its support is not finite. Such priors do not
necessarily result in a proper posterior distribution. Northrop and
Attalides (2016) consider the issue of posterior propriety in Bayesian
extreme value analyses. In most of improper priors below the prior for
the scale parameter \\\sigma\\ is taken to be \\1/\sigma\\, i.e. a flat
prior for \\\log \sigma\\. Here we denote the scale parameter of the GP
distribution by \\\sigma\\, whereas we use \\\sigma_u\\ in the revdbayes
vignette.

For all in-built priors the arguments `min_xi` and `max_xi` may be
supplied by the user. The prior density is set to zero for any value of
the shape parameter \\\xi\\ that is outside (`min_xi`, `max_xi`). This
will override the default values of `min_xi` and `max_xi` in the named
priors detailed above.

**Extreme value priors.** It is typical to use either `prior = "prob"`
([`gev_prob`](https://github.com/paulnorthrop/revdbayes/reference/gev_prob.md))
or `prior = "quant"`
([`gev_quant`](https://github.com/paulnorthrop/revdbayes/reference/gev_quant.md))
to set an informative prior and one of the other prior (or a
user-supplied function) otherwise. The names of the in-built extreme
value priors set using `prior` and details of hyperparameters are:

- `"prob"`. A prior for GEV parameters \\(\mu, \sigma, \xi)\\ based on
  Crowder (1992). See
  [`gev_prob`](https://github.com/paulnorthrop/revdbayes/reference/gev_prob.md)
  for details. See also Northrop et al. (2017) and Stephenson (2016).

- `"quant"`. A prior for GEV parameters \\(\mu, \sigma, \xi)\\ based on
  Coles and Tawn (1996). See
  [`gev_quant`](https://github.com/paulnorthrop/revdbayes/reference/gev_quant.md)
  for details.

- `"norm"`.

  For `model = "gp"`: (\\\log \sigma, \xi\\), is bivariate normal with
  mean `mean` (a numeric vector of length 2) and covariance matrix `cov`
  (a symmetric positive definite 2 by 2 matrix).

  For `model = "gev"`: (\\\mu, \log \sigma, \xi\\), is trivariate normal
  with mean `mean` (a numeric vector of length 3) and covariance matrix
  `cov` (a symmetric positive definite 3 by 3 matrix).

- `"loglognorm"`. For `model = "gev"` only: (\\\log \mu, \log \sigma,
  \xi\\), is trivariate normal with mean `mean` (a numeric vector of
  length 3) and covariance matrix `cov` (a symmetric positive definite 3
  by 3 matrix).

- `"mdi"`.

  For `model = "gp"`: (an extended version of) the maximal data
  information (MDI) prior, that is, \$\$\pi(\sigma, \xi) = \sigma^{-1}
  \exp\[-a(\xi + 1)\], {\rm ~for~} \sigma \> 0, \xi \geq -1, a \geq
  0.\$\$ The value of \\a\\ is set using the argument `a`. The default
  value is \\a = 1\\, which gives the MDI prior.

  For `model = "gev"`: (an extended version of) the maximal data
  information (MDI) prior, that is, \$\$\pi(\mu, \sigma, \xi) =
  \sigma^{-1} \exp\[-a(\xi + 1)\], {\rm ~for~} \sigma \> 0, \xi \geq -1,
  a \geq 0.\$\$ The value of \\a\\ is set using the argument `a`. The
  default value is \\a = \gamma\\, where \\\gamma = 0.57721\\ is Euler's
  constant, which gives the MDI prior.

  For each of these cases \\\xi\\ must be is bounded below *a priori*
  for the posterior to be proper (Northrop and Attalides, 2016). An
  argument for the bound \\\xi \geq -1\\ is that for \\\xi \< -1\\ the
  GP (GEV) likelihood is unbounded above as \\-\sigma/\xi\\ (\\\mu -
  \sigma/\xi\\)) approaches the sample maximum. In maximum likelihood
  estimation of GP parameters (Grimshaw, 1993) and GEV parameters a
  local maximum of the likelihood is sought on the region \\\sigma \> 0,
  \xi \geq -1\\.

- `"flat"`.

  For `model = "gp"`: a flat prior for \\\xi\\ (and for \\\log
  \sigma\\): \$\$\pi(\sigma, \xi) = \sigma^{-1}, {\rm ~for~} \sigma \>
  0.\$\$

  For `model = "gev"`: a flat prior for \\\xi\\ (and for \\\mu\\ and
  \\\log \sigma\\): \$\$\pi(\mu, \sigma, \xi) = \sigma^{-1}, {\rm ~for~}
  \sigma \> 0.\$\$

  These are the conditional right Haar priors of Jewson et al. (2025).

- `"flatflat"`.

  For `model = "gp"`: flat priors for \\\sigma\\ and \\\xi\\:
  \$\$\pi(\sigma, \xi) = 1, {\rm ~for~} \sigma \> 0.\$\$

  For `model = "gev"`: flat priors for \\\mu\\, \\\sigma\\ and \\\xi\\:
  \$\$\pi(\mu, \sigma, \xi) = 1, {\rm ~for~} \sigma \> 0.\$\$

  Therefore, the posterior is proportional to the likelihood.

- `"jeffreys"`. For `model = "gp"` only: the Jeffreys prior (Castellanos
  and Cabras, 2007): \$\$\pi(\sigma, \xi) =
  \sigma^{-1}(1+\xi)^{-1}(1+2\xi)^{-1/2}, {\rm ~for~} \sigma \> 0, \xi
  \> -1 / 2.\$\$

  In the GEV case the Jeffreys prior doesn't yield a proper posterior
  for any sample size. See Northrop and Attalides (2016) for details.

- `"beta"`. For `model = "gp"`: a beta-type(p, q) prior is used for xi
  on the interval (`min_xi`, `max_xi`): \$\$\pi(\sigma, \xi) =
  \sigma^{-1} (\xi - {\min}\_{\xi}) ^ {p-1} ({\max}\_{\xi} - \xi) ^
  {q-1}, {\rm ~for~} {\min}\_{\xi} \< \xi \< {\max}\_{\xi}.\$\$

  For `model = "gev"`: similarly ... \$\$\pi(\mu, \sigma, \xi) =
  \sigma^{-1} (\xi - {\min}\_{\xi}) ^ {p-1} ({\max}\_{\xi} - \xi) ^
  {q-1}, {\rm ~for~} {\min}\_{\xi} \< \xi \< {\max}\_{\xi}.\$\$

  The argument `pq` is a vector containing `c(p,q)`. The default
  settings for this prior are `p = 6, q = 9` and
  `min_xi = -1/2, max_xi = 1/2`, which corresponds to the prior for
  \\\xi\\ proposed in Martins and Stedinger (2000, 2001).

## References

Castellanos, E. M. and Cabras, S. (2007) A default Bayesian procedure
for the generalized Pareto distribution. *Journal of Statistical
Planning and Inference* **137(2)**, 473-483.
[doi:10.1016/j.jspi.2006.01.006](https://doi.org/10.1016/j.jspi.2006.01.006)
.

Coles, S. G. and Tawn, J. A. (1996) A Bayesian analysis of extreme
rainfall data. *Appl. Statist.*, **45**, 463-478.

Crowder, M. (1992) Bayesian priors based on parameter transformation
using the distribution function *Ann. Inst. Statist. Math.*, **44**,
405-416. <https://link.springer.com/article/10.1007/BF00050695>.

Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates for the
Generalized Pareto Distribution. *Technometrics*, **35(2)**, 185-191.
[doi:10.1080/00401706.1993.10485040](https://doi.org/10.1080/00401706.1993.10485040)
.

Hosking, J. R. M. and Wallis, J. R. (1987) Parameter and Quantile
Estimation for the Generalized Pareto Distribution. *Technometrics*,
**29(3)**, 339-349.
[doi:10.2307/1269343](https://doi.org/10.2307/1269343) .

Jewson, S., Sweeting, T., and Jewson, L. (2025) Reducing reliability
bias in assessments of extreme weather risk using calibrating priors,
*Adv. Stat. Clim. Meteorol. Oceanogr.*, **11**, 1-22,
[doi:10.5194/ascmo-11-1-2025](https://doi.org/10.5194/ascmo-11-1-2025) .

Martins, E. S. and Stedinger, J. R. (2000) Generalized maximum
likelihood generalized extreme value quantile estimators for hydrologic
data, *Water Resources Research*, **36(3)**, 737-744.
[doi:10.1029/1999WR900330](https://doi.org/10.1029/1999WR900330) .

Martins, E. S. and Stedinger, J. R. (2001) Generalized maximum
likelihood Pareto-Poisson estimators for partial duration series, *Water
Resources Research*, **37(10)**, 2551-2557.
[doi:10.1029/2001WR000367](https://doi.org/10.1029/2001WR000367) .

Northrop, P.J. and Attalides, N. (2016) Posterior propriety in Bayesian
extreme value analyses using reference priors *Statistica Sinica*,
**26**(2), 721–743
[doi:10.5705/ss.2014.034](https://doi.org/10.5705/ss.2014.034) .

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

[`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
and
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md)
for sampling from an extreme value posterior distribution.

[`create_prior_xptr`](https://github.com/paulnorthrop/revdbayes/reference/create_prior_xptr.md)
for creating an external pointer to a C++ function to evaluate the
log-prior density.

[`rprior_prob`](https://github.com/paulnorthrop/revdbayes/reference/rprior_prob.md)
and
[`rprior_quant`](https://github.com/paulnorthrop/revdbayes/reference/rprior_quant.md)
for sampling from informative prior distributions for GEV parameters.

[`gp_norm`](https://github.com/paulnorthrop/revdbayes/reference/gp_norm.md),
[`gp_mdi`](https://github.com/paulnorthrop/revdbayes/reference/gp_mdi.md),
[`gp_flat`](https://github.com/paulnorthrop/revdbayes/reference/gp_flat.md),
[`gp_flatflat`](https://github.com/paulnorthrop/revdbayes/reference/gp_flatflat.md),
[`gp_jeffreys`](https://github.com/paulnorthrop/revdbayes/reference/gp_jeffreys.md),
[`gp_beta`](https://github.com/paulnorthrop/revdbayes/reference/gp_beta.md)
to see the arguments for priors for GP parameters.

[`gev_norm`](https://github.com/paulnorthrop/revdbayes/reference/gev_norm.md),
[`gev_loglognorm`](https://github.com/paulnorthrop/revdbayes/reference/gev_loglognorm.md),
[`gev_mdi`](https://github.com/paulnorthrop/revdbayes/reference/gev_mdi.md),
[`gev_flat`](https://github.com/paulnorthrop/revdbayes/reference/gev_flat.md),
[`gev_flatflat`](https://github.com/paulnorthrop/revdbayes/reference/gev_flatflat.md),
[`gev_beta`](https://github.com/paulnorthrop/revdbayes/reference/gev_beta.md),
[`gev_prob`](https://github.com/paulnorthrop/revdbayes/reference/gev_prob.md),
[`gev_quant`](https://github.com/paulnorthrop/revdbayes/reference/gev_quant.md)
to see the arguments for priors for GEV parameters.

## Examples

``` r
# Normal prior for GEV parameters (mu, log(sigma), xi).
mat <- diag(c(10000, 10000, 100))
pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
pn
#> $prior
#> [1] "gev_norm"
#> 
#> $mean
#> [1] 0 0 0
#> 
#> $min_xi
#> [1] -Inf
#> 
#> $max_xi
#> [1] Inf
#> 
#> $icov
#> [1] 1e-04 0e+00 0e+00 1e-04 0e+00 1e-02
#> 
#> $trendsd
#> [1] 0
#> 
#> attr(,"class")
#> [1] "evprior"
#> attr(,"model")
#> [1] "gev"

# Prior for GP parameters with flat prior for xi on (-1, infinity).
fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
fp
#> $prior
#> [1] "gp_flat"
#> 
#> $min_xi
#> [1] -1
#> 
#> $max_xi
#> [1] Inf
#> 
#> $trendsd
#> [1] 0
#> 
#> attr(,"class")
#> [1] "evprior"
#> attr(,"model")
#> [1] "gp"

# A user-defined prior (see the vignette for details).
u_prior_fn <- function(x, ab) {
  if (x[1] <= 0 | x[2] <= -1 | x[2] >= 1) {
    return(-Inf)
  }
  return(-log(x[1]) + (ab[1] - 1) * log(1 + x[2]) +
         (ab[2] - 1) * log(1 - x[2]))
}
up <- set_prior(prior = u_prior_fn, ab = c(2, 2), model = "gp")

# A user-defined prior using a pointer to a C++ function
ptr_gp_flat <- create_prior_xptr("gp_flat")
u_prior_ptr <- set_prior(prior = ptr_gp_flat, model = "gp")
```
