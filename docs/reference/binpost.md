# Random sampling from a binomial posterior distribution

Samples from the posterior distribution of the probability \\p\\ of a
binomial distribution.

## Usage

``` r
binpost(n, prior, ds_bin, param = c("logit", "p"))
```

## Arguments

- n:

  A numeric scalar. The size of posterior sample required.

- prior:

  A function to evaluate the prior, created by
  [`set_bin_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_bin_prior.md).

- ds_bin:

  A numeric list. Sufficient statistics for inference about a binomial
  probability \\p\\. Contains

  - `n_raw` : number of raw observations.

  - `m` : number of threshold exceedances.

- param:

  A character scalar. Only relevant if `prior$prior` is a
  (user-supplied) R function. `param` specifies the parameterization of
  the posterior distribution that
  [`ru`](https://github.com/paulnorthrop/rust/reference/ru.html) uses
  for sampling.

  If `param = "p"` the original parameterization \\p\\ is used.

  If `param = "logit"` (the default) then
  [`ru`](https://github.com/paulnorthrop/rust/reference/ru.html) samples
  from the posterior for the logit of \\p\\, before transforming back to
  the \\p\\-scale.

  The latter tends to make the optimizations involved in the
  ratio-of-uniforms algorithm more stable and to increase the
  probability of acceptance, but at the expense of slower function
  evaluations.

## Value

An object (list) of class `"binpost"` with components

- `bin_sim_vals`::

  An `n` by 1 numeric matrix of values simulated from the posterior for
  the binomial probability \\p\\

- `bin_logf`::

  A function returning the log-posterior for \\p\\.

- `bin_logf_args`::

  A list of arguments to `bin_logf`.

If `prior$prior` is a (user-supplied) R function then this list also
contains `ru_object` the object of class `"ru"` returned by
[`ru`](https://github.com/paulnorthrop/rust/reference/ru.html).

## Details

If `prior$prior == "bin_beta"` then the posterior for \\p\\ is a beta
distribution so [`rbeta`](https://rdrr.io/r/stats/Beta.html) is used to
sample from the posterior.

If `prior$prior == "bin_mdi"` then rejection sampling is used to sample
from the posterior with an envelope function equal to the density of a
beta(`ds$m` + 1, `ds$n_raw - ds$m` + 1) density.

If `prior$prior == "bin_northrop"` then rejection sampling is used to
sample from the posterior with an envelope function equal to the
posterior density that results from using a Haldane prior.

If `prior$prior` is a (user-supplied) R function then
[`ru`](https://github.com/paulnorthrop/rust/reference/ru.html) is used
to sample from the posterior using the generalised ratio-of-uniforms
method.

## See also

[`set_bin_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_bin_prior.md)
for setting a prior distribution for the binomial probability \\p\\.

## Examples

``` r
u <- quantile(gom, probs = 0.65)
ds_bin <- list()
ds_bin$n_raw <- length(gom)
ds_bin$m <- sum(gom > u)
bp <- set_bin_prior(prior = "jeffreys")
temp <- binpost(n = 1000, prior = bp, ds_bin = ds_bin)
graphics::hist(temp$bin_sim_vals, prob = TRUE)


# Setting a beta prior (Jeffreys in this case) by hand
beta_prior_fn <- function(p, ab) {
  return(stats::dbeta(p, shape1 = ab[1], shape2 = ab[2], log = TRUE))
}
jeffreys <- set_bin_prior(beta_prior_fn, ab = c(1 / 2, 1 / 2))
temp <- binpost(n = 1000, prior = jeffreys, ds_bin = ds_bin)
```
