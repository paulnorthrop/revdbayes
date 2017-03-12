# revdbayes 1.1.0

## New features

* A new vignette (Posterior Predictive Extreme Value Inference using the revdbayes   
  Package) provides an overview of most of the new features. Run 
  browseVignettes("revdbayes") to access.

* S3 `predict()` method for class 'evpost' performs predictive inference 
  about the largest observation observed in N years, returning an object
  of class `evpred`.
  
* S3 `plot()` for the `evpred` object returned by `predict.evpost`.

* S3 `pp_check()` method for class 'evpost' performs posterior predictive 
  checks using the bayesplot package.

* Interface to the bayesplot package added in the S3 `plot.evpost` method.

* `model = bingp` can now be supplied to `rpost()` to add inferences about the
  probability of threshold exceedance to inferences about threshold excesses 
  based on the Generalised Pareto (GP) model.  `set_bin_prior()` can be used to 
  set a prior for this probability.

* `rprior_quant()`: to simulate from the prior distribution for GEV parameters 
  proposed in Coles and Tawn (1996) [A Bayesian analysis of extreme rainfall 
  data. Appl. Statist., 45, 463-478], based on independent gamma priors for 
  differences between quantiles.  
   
* `prior_prob()`: to simulate from the prior distribution for GEV parameters
  based on Crowder (1992), in which independent beta priors are specified for 
  ratios of probabilities (which is equivalent to a Dirichlet prior on 
  differences between these probabilities).

## Bug fixes and minor improvements

* The spurious warning messages relating to checking that the model argument
  to `rpost()` is consistent with the prior set using `set-prior()` have been 
  corrected.  These occurred when `model = "pp"` or `model = "os"`.
  
* The hyperparameter in the MDI prior was `a` in the documentation and `a_mdi`
  in the code.  Now it is `a` everywhere.
  
* In `set_prior` with `prior = "beta"` parameter vector `ab` has been 
  corrected to `pq`.
  
* In the documentation of `rpost()` the description of the argument `noy` 
  has been corrected.
  
* Package spatstat removed from the Imports field in description to avoid
  NOTE in CRAN checks.  

