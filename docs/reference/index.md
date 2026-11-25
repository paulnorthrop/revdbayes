# Package index

## Package

- [`revdbayes`](https://github.com/paulnorthrop/revdbayes/reference/revdbayes-package.md)
  [`revdbayes-package`](https://github.com/paulnorthrop/revdbayes/reference/revdbayes-package.md)
  : revdbayes: Ratio-of-Uniforms Sampling for Bayesian Extreme Value
  Analysis

## Priors

- [`set_prior()`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md)
  : Construction of prior distributions for extreme value model
  parameters
- [`set_bin_prior()`](https://github.com/paulnorthrop/revdbayes/reference/set_bin_prior.md)
  : Construction of a prior distribution for a binomial probability
  \\p\\
- [`gev_beta()`](https://github.com/paulnorthrop/revdbayes/reference/gev_beta.md)
  : Beta-type prior for GEV shape parameter \\\xi\\
- [`gev_flat()`](https://github.com/paulnorthrop/revdbayes/reference/gev_flat.md)
  : Flat prior for GEV parameters (\\\mu, log \sigma, \xi\\)
- [`gev_flatflat()`](https://github.com/paulnorthrop/revdbayes/reference/gev_flatflat.md)
  : Flat prior for GEV parameters (\\\mu, \sigma, \xi\\)
- [`gev_loglognorm()`](https://github.com/paulnorthrop/revdbayes/reference/gev_loglognorm.md)
  : Trivariate normal prior for GEV parameters (\\log \mu, log \sigma,
  \xi\\)
- [`gev_mdi()`](https://github.com/paulnorthrop/revdbayes/reference/gev_mdi.md)
  : Maximal data information (MDI) prior for GEV parameters (\\\mu,
  \sigma, \xi\\)
- [`gev_norm()`](https://github.com/paulnorthrop/revdbayes/reference/gev_norm.md)
  : Trivariate normal prior for GEV parameters (\\\mu, log \sigma,
  \xi\\)
- [`gev_prob()`](https://github.com/paulnorthrop/revdbayes/reference/gev_prob.md)
  : Informative GEV prior on a probability scale
- [`gev_quant()`](https://github.com/paulnorthrop/revdbayes/reference/gev_quant.md)
  : Informative GEV prior on a quantile scale
- [`gp_beta()`](https://github.com/paulnorthrop/revdbayes/reference/gp_beta.md)
  : Beta-type prior for GP shape parameter \\\xi\\
- [`gp_flat()`](https://github.com/paulnorthrop/revdbayes/reference/gp_flat.md)
  : Flat prior for GP parameters (\\log \sigma, \xi\\)
- [`gp_flatflat()`](https://github.com/paulnorthrop/revdbayes/reference/gp_flatflat.md)
  : Flat prior for GP parameters (\\\sigma, \xi\\)
- [`gp_jeffreys()`](https://github.com/paulnorthrop/revdbayes/reference/gp_jeffreys.md)
  : Jeffreys prior for GP parameters (\\\sigma, \xi\\)
- [`gp_mdi()`](https://github.com/paulnorthrop/revdbayes/reference/gp_mdi.md)
  : Maximal data information (MDI) prior for GP parameters (\\\sigma,
  \xi\\)
- [`gp_norm()`](https://github.com/paulnorthrop/revdbayes/reference/gp_norm.md)
  : Bivariate normal prior for GP parameters (\\log \sigma, \xi\\)
- [`rprior_prob()`](https://github.com/paulnorthrop/revdbayes/reference/rprior_prob.md)
  : Prior simulation of GEV parameters - prior on probability scale
- [`rprior_quant()`](https://github.com/paulnorthrop/revdbayes/reference/rprior_quant.md)
  : Prior simulation of GEV parameters - prior on quantile scale
- [`rDir()`](https://github.com/paulnorthrop/revdbayes/reference/rDir.md)
  : Simulation from a Dirichlet distribution

## Posterior sampling

- [`rpost()`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
  : Random sampling from extreme value posterior distributions
- [`rpost_rcpp()`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md)
  : Random sampling from extreme value posterior distributions
- [`kgaps_post()`](https://github.com/paulnorthrop/revdbayes/reference/kgaps_post.md)
  : Random sampling from K-gaps posterior distribution
- [`dgaps_post()`](https://github.com/paulnorthrop/revdbayes/reference/dgaps_post.md)
  : Random sampling from D-gaps posterior distribution
- [`binpost()`](https://github.com/paulnorthrop/revdbayes/reference/binpost.md)
  : Random sampling from a binomial posterior distribution
- [`wbinpost()`](https://github.com/paulnorthrop/revdbayes/reference/wbinpost.md)
  : Random sampling from a binomial posterior distribution, using
  weights

## Predictive inference and checking

- [`predict(`*`<evpost>`*`)`](https://github.com/paulnorthrop/revdbayes/reference/predict.evpost.md)
  : Predictive inference for the largest value observed in \\N\\ years.
- [`pp_check(`*`<evpost>`*`)`](https://github.com/paulnorthrop/revdbayes/reference/pp_check.evpost.md)
  : Posterior predictive checks for an evpost object

## Distributions

- [`dgev()`](https://github.com/paulnorthrop/revdbayes/reference/gev.md)
  [`pgev()`](https://github.com/paulnorthrop/revdbayes/reference/gev.md)
  [`qgev()`](https://github.com/paulnorthrop/revdbayes/reference/gev.md)
  [`rgev()`](https://github.com/paulnorthrop/revdbayes/reference/gev.md)
  : The Generalised Extreme Value Distribution
- [`dgp()`](https://github.com/paulnorthrop/revdbayes/reference/gp.md)
  [`pgp()`](https://github.com/paulnorthrop/revdbayes/reference/gp.md)
  [`qgp()`](https://github.com/paulnorthrop/revdbayes/reference/gp.md)
  [`rgp()`](https://github.com/paulnorthrop/revdbayes/reference/gp.md) :
  The Generalised Pareto Distribution

## Example data

- [`gom`](https://github.com/paulnorthrop/revdbayes/reference/gom.md) :
  Storm peak significant wave heights from the Gulf of Mexico
- [`newlyn`](https://github.com/paulnorthrop/revdbayes/reference/newlyn.md)
  : Newlyn sea surges
- [`oxford`](https://github.com/paulnorthrop/revdbayes/reference/oxford.md)
  : Annual Maximum Temperatures at Oxford
- [`portpirie`](https://github.com/paulnorthrop/revdbayes/reference/portpirie.md)
  : Annual Maximum Sea Levels at Port Pirie, South Australia
- [`rainfall`](https://github.com/paulnorthrop/revdbayes/reference/rainfall.md)
  : Daily Aggregate Rainfall
- [`venice`](https://github.com/paulnorthrop/revdbayes/reference/venice.md)
  : Largest Sea Levels in Venice

## Generics

- [`plot(`*`<evpost>`*`)`](https://github.com/paulnorthrop/revdbayes/reference/plot.evpost.md)
  : Plot diagnostics for an evpost object
- [`plot(`*`<evpred>`*`)`](https://github.com/paulnorthrop/revdbayes/reference/plot.evpred.md)
  : Plot diagnostics for an evpred object
- [`print(`*`<evpost>`*`)`](https://github.com/paulnorthrop/revdbayes/reference/print.evpost.md)
  : Print method for objects of class "evpost"
- [`summary(`*`<evpost>`*`)`](https://github.com/paulnorthrop/revdbayes/reference/summary.evpost.md)
  : Summarizing an evpost object
- [`print(`*`<summary.evpost>`*`)`](https://github.com/paulnorthrop/revdbayes/reference/print.summary.evpost.md)
  : Print method for objects of class "summary.evpost"

## Miscellaneous

- [`quantile_to_gev()`](https://github.com/paulnorthrop/revdbayes/reference/quantile_to_gev.md)
  : Converts quantiles to GEV parameters
- [`create_prior_xptr()`](https://github.com/paulnorthrop/revdbayes/reference/create_prior_xptr.md)
  : Create an external pointer to a C++ prior
- [`grimshaw_gp_mle()`](https://github.com/paulnorthrop/revdbayes/reference/grimshaw_gp_mle.md)
  : Maximum likelihood estimation of generalised Pareto parameters
- [`gp_pwm()`](https://github.com/paulnorthrop/revdbayes/reference/gp_pwm.md)
  : Probability-weighted moments estimation of generalised Pareto
  parameters
- [`gp_lrs()`](https://github.com/paulnorthrop/revdbayes/reference/gp_lrs.md)
  : Linear Combinations of Ratios of Spacings estimation of generalised
  Pareto parameters
