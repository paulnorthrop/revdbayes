# Random sampling from K-gaps posterior distribution

Uses the
[`rust`](https://github.com/paulnorthrop/rust/reference/rust-package.html)
package to simulate from the posterior distribution of the extremal
index \\\theta\\ based on the K-gaps model for threshold interexceedance
times of Suveges and Davison (2010).

## Usage

``` r
kgaps_post(
  data,
  thresh,
  k = 1,
  n = 1000,
  inc_cens = TRUE,
  alpha = 1,
  beta = 1,
  param = c("logit", "theta"),
  use_rcpp = TRUE
)
```

## Arguments

- data:

  A numeric vector or numeric matrix of raw data. If `data` is a matrix
  then the log-likelihood is constructed as the sum of (independent)
  contributions from different columns. A common situation is where each
  column relates to a different year.

  If `data` contains missing values then
  [`split_by_NAs`](https://github.com/paulnorthrop/exdex/reference/split_by_NAs.html)
  is used to divide the data further into sequences of non-missing
  values, stored in different columns in a matrix. Again, the
  log-likelihood is constructed as a sum of contributions from different
  columns.

- thresh:

  A numeric scalar. Extreme value threshold applied to data.

- k:

  A numeric scalar. Run parameter \\K\\, as defined in Suveges and
  Davison (2010). Threshold inter-exceedances times that are not larger
  than `k` units are assigned to the same cluster, resulting in a
  \\K\\-gap equal to zero. Specifically, the \\K\\-gap \\S\\
  corresponding to an inter-exceedance time of \\T\\ is given by \\S =
  \max(T - K, 0)\\.

- n:

  A numeric scalar. The size of posterior sample required.

- inc_cens:

  A logical scalar indicating whether or not to include contributions
  from right-censored inter-exceedance times, relating to the first and
  last observations. It is known that these times are greater than or
  equal to the time observed. If `data` has multiple columns then there
  will be right-censored first and last inter-exceedance times for each
  column. See also the **Details** section of
  [`kgaps`](https://github.com/paulnorthrop/exdex/reference/kgaps.html).

- alpha, beta:

  Positive numeric scalars. Parameters of a beta(\\\alpha\\, \\\beta\\)
  prior for \\\theta\\.

- param:

  A character scalar. If `param = "logit"` (the default) then we
  simulate from the posterior distribution of \\\phi = \log(\theta /
  (1-\theta))\\ and then transform back to the \\\theta\\-scale. If
  `param = "theta"` then we simulate directly from the posterior
  distribution of \\\theta\\, unless the sample K-gaps are all equal to
  zero or all positive, when we revert to `param = "logit"`. This is to
  avoid sampling directly from a posterior with mode equal to 0 or 1.

- use_rcpp:

  A logical scalar. If `TRUE` (the default) the rust function
  [`ru_rcpp`](https://github.com/paulnorthrop/rust/reference/ru_rcpp.html)
  is used for posterior simulation. If `FALSE` the (slower) function
  [`ru`](https://github.com/paulnorthrop/rust/reference/ru.html) is
  used.

## Value

An object (list) of class `"evpost"`, which has the same structure as an
object of class `"ru"` returned from
[`ru`](https://github.com/paulnorthrop/rust/reference/ru.html). In
addition this list contains

- `call`: The call to `kgaps()`.

- `model`: The character scalar `"kgaps"`.

- `thresh`: The argument `thresh`.

- `ss`: The sufficient statistics for the K-gaps likelihood, as
  calculated by
  [`kgaps_stat`](https://github.com/paulnorthrop/exdex/reference/kgaps_stat.html).

## Details

A beta(\\\alpha\\, \\\beta\\) prior distribution is used for \\\theta\\
so that the posterior from which values are simulated is proportional to
\$\$\theta ^ {2 N_1 + \alpha - 1} (1 - \theta) ^ {N_0 + \beta - 1}
\exp\\- \theta q (S_0 + \cdots + S_N)\\.\$\$ See
[`kgaps_stat`](https://github.com/paulnorthrop/exdex/reference/kgaps_stat.html)
for a description of the variables involved in the contribution of the
likelihood to this expression.

The [`ru`](https://github.com/paulnorthrop/rust/reference/ru.html)
function in the
[`rust`](https://github.com/paulnorthrop/rust/reference/rust-package.html)
package simulates from this posterior distribution using the generalised
ratio-of-uniforms distribution. To improve the probability of
acceptance, and to ensure that the simulation will work even in extreme
cases where the posterior density of \\\theta\\ is unbounded as
\\\theta\\ approaches 0 or 1, we simulate from the posterior
distribution of \\\phi = \log(\theta / (1-\theta))\\ and then transform
back to the \\\theta\\-scale.

## References

Suveges, M. and Davison, A. C. (2010) Model misspecification in peaks
over threshold analysis, *The Annals of Applied Statistics*, **4**(1),
203-221. [doi:10.1214/09-AOAS292](https://doi.org/10.1214/09-AOAS292)

## See also

[`ru`](https://github.com/paulnorthrop/rust/reference/ru.html) for the
form of the object returned by `kgaps_post`.

[`dgaps_post`](https://github.com/paulnorthrop/revdbayes/reference/dgaps_post.md)
for Bayesian inference about the extremal index \\\theta\\ using the
\\D\\-gaps model.

## Examples

``` r
# Newlyn sea surges

thresh <- quantile(newlyn, probs = 0.90)
k_postsim <- kgaps_post(newlyn, thresh)
plot(k_postsim)


### Cheeseboro wind gusts

k_postsim <- kgaps_post(exdex::cheeseboro, thresh = 45, k = 3)
plot(k_postsim)
```
