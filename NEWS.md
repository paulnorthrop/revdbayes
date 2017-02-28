# revdbayes 1.0.0.9000

## New features

* `model = bingp` can now be supplied to `rpost()` to add inferences about the
  probabiliy of threshold exceedance to inferences about threshold excesses 
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
  
* Functions `dpred()`, `ppred()`, `qpred()` and `rpred()` to perform predictive
  inference about the largest observation observed in N years.

## Bug fixes and minor improvements

* The spurious warning messages relating to checking that the model argument
  to `rpost()` is consistent with the prior set using `set-prior()` have been 
  corrected.  These occurred when `model = "pp"` or `model = "os"`.
  
* The hyperparameter in the MDI prior was `a` in the documentation and `a_mdi`
  in the code.  Now it is `a` everywhere.
  
* In the documentation of `rpost()` the description of the argument `noy` 
  has been corrected.

