# Random sampling from a binomial posterior distribution, using weights

Samples from the posterior distribution of the probability \\p\\ of a
binomial distribution. User-supplied weights are applied to each
observation when constructing the log-likelihood.

## Usage

``` r
wbinpost(n, prior, ds_bin)
```

## Arguments

- n:

  A numeric scalar. The size of posterior sample required.

- prior:

  A function to evaluate the prior, created by
  [`set_bin_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_bin_prior.md).
  `prior$prior` must be `"bin_beta"`.

- ds_bin:

  A numeric list. Sufficient statistics for inference about the binomial
  probability \\p\\. Contains

  - `sf` : a logical vector of success (`TRUE`) and failure (`FALSE`)
    indicators.

  - `w` : a numeric vector of length `length(sf)` containing the values
    by which to multiply the observations when constructing the
    log-likelihood.

## Value

An object (list) of class `"binpost"` with components

- `bin_sim_vals`::

  An `n` by 1 numeric matrix of values simulated from the posterior for
  the binomial probability \\p\\

- `bin_logf`::

  A function returning the log-posterior for \\p\\.

- `bin_logf_args`::

  A list of arguments to `bin_logf`.

## Details

For `prior$prior == "bin_beta"` the posterior for \\p\\ is a beta
distribution so [`rbeta`](https://rdrr.io/r/stats/Beta.html) is used to
sample from the posterior.

## See also

[`set_bin_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_bin_prior.md)
for setting a prior distribution for the binomial probability \\p\\.

## Examples

``` r
u <- quantile(gom, probs = 0.65)
ds_bin <- list(sf = gom > u, w = rep(1, length(gom)))
bp <- set_bin_prior(prior = "jeffreys")
temp <- wbinpost(n = 1000, prior = bp, ds_bin = ds_bin)
graphics::hist(temp$bin_sim_vals, prob = TRUE)
```
