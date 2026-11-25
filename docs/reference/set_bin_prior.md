# Construction of a prior distribution for a binomial probability \\p\\

Constructs a prior distribution for use as the argument `bin_prior` in
[`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
or in
[`binpost`](https://github.com/paulnorthrop/revdbayes/reference/binpost.md).
The user can choose from a list of in-built priors or specify their own
prior function, returning the **log** of the prior density, using an R
function and arguments for hyperparameters.

## Usage

``` r
set_bin_prior(
  prior = c("jeffreys", "laplace", "haldane", "beta", "mdi", "northrop"),
  ...
)
```

## Arguments

- prior:

  Either

  - An R function that returns the value of the log of the prior density
    (see **Examples**), or

  - A character string giving the name of the prior for \\p\\. See
    **Details** for a list of priors available.

- ...:

  Further arguments to be passed to the user-supplied or in-built prior
  function. For the latter this is only relevant if `prior = "beta"`,
  when `ab` can be passed. See **Details**.

## Value

A list of class `"binprior"`. The first component is the name of the
input prior. Apart from the MDI prior this will be "beta", in which case
the other component of the list is a vector of length two giving the
corresponding values of the beta parameters.

## Details

**Binomial priors.** The names of the binomial priors set using
`bin_prior` are:

- `"jeffreys"`: the *Jeffreys* beta(1/2, 1/2) prior.

- `"laplace"`: the *Bayes-Laplace* beta(1, 1) prior.

- `"haldane"`: the *Haldane* beta(0, 0) prior.

- `"beta"`: a beta(\\\alpha, \beta\\) prior. The argument `ab` is a
  vector containing `c`(\\\alpha, \beta\\). The default is
  `ab = c(1, 1)`.

- `"mdi"`: the MDI prior \\\pi(p) = 1.6186 p^p (1-p)^{1-p}\\, for \\0 \<
  p \< 1.\\

- `"northrop"`: the improper prior
  \\\pi(p)=\\-\ln(1-p)\\^{-1}(1-p)^{-1}\\, for \\0 \< p \< 1.\\

Apart from the last two priors these are all beta distributions.

## See also

[`binpost`](https://github.com/paulnorthrop/revdbayes/reference/binpost.md)
for sampling from a binomial posterior distribution.

## Examples

``` r
bp <- set_bin_prior(prior = "jeffreys")

# Setting the Jeffreys prior by hand
beta_prior_fn <- function(p, ab) {
  return(stats::dbeta(p, shape1 = ab[1], shape2 = ab[2], log = TRUE))
}
jeffreys <- set_bin_prior(beta_prior_fn, ab = c(1 / 2, 1 / 2))
```
