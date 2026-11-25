# revdbayes: Ratio-of-Uniforms Sampling for Bayesian Extreme Value Analysis

Uses the multivariate generalized ratio-of-uniforms method to simulate
random samples from the posterior distributions commonly encountered in
Bayesian extreme value analyses.

## Details

The main functions in the revdbayes package are
[`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
and
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md),
which simulate random samples from the posterior distribution of extreme
value model parameters using the functions
[`ru`](https://github.com/paulnorthrop/rust/reference/ru.html) and
[`ru_rcpp`](https://github.com/paulnorthrop/rust/reference/ru_rcpp.html)
from the rust package, respectively. The user chooses the extreme value
model, the prior density for the parameters and provides the data. There
are options to improve the probability of acceptance of the
ratio-of-uniforms algorithm by working with transformation of the model
parameters.

The functions
[`kgaps_post`](https://github.com/paulnorthrop/revdbayes/reference/kgaps_post.md)
and
[`dgaps_post`](https://github.com/paulnorthrop/revdbayes/reference/dgaps_post.md)
simulate from the posterior distribution of the extremal index
\\\theta\\ based on the K-gaps model for threshold interexceedance times
of Suveges and Davison (2010) and the similar D-gaps model of Holesovsky
and Fusek (2020). See also Attalides (2015).

See
[`vignette("revdbayes-a-vignette", package = "revdbayes")`](https://github.com/paulnorthrop/revdbayes/articles/revdbayes-a-vignette.md)
for an overview of the package and
[`vignette("revdbayes-b-using-rcpp-vignette", package = "revdbayes")`](https://github.com/paulnorthrop/revdbayes/articles/revdbayes-b-using-rcpp-vignette.md)
for an illustration of the improvements in efficiency produced using the
Rcpp package. See
[`vignette("revdbayes-c-predictive-vignette", package = "revdbayes")`](https://github.com/paulnorthrop/revdbayes/articles/revdbayes-c-predictive-vignette.md)
for an outline of how to use revdbayes to perform posterior predictive
extreme value inference and
[`vignette("revdbayes-d-kgaps-vignette", package = "revdbayes")`](https://github.com/paulnorthrop/revdbayes/articles/revdbayes-d-kgaps-vignette.md)
considers Bayesian inference for the extremal index \\\theta\\ using
threshold inter-exceedance times.

## References

Holesovsky, J. and Fusek, M. Estimation of the extremal index using
censored distributions. Extremes 23, 197-213 (2020).
[doi:10.1007/s10687-020-00374-3](https://doi.org/10.1007/s10687-020-00374-3)

Northrop, P. J. (2016). rust: Ratio-of-Uniforms Simulation with
Transformation. R package version 1.2.2.
<https://cran.r-project.org/package=rust>.

Suveges, M. and Davison, A. C. (2010) Model misspecification in peaks
over threshold analysis, *The Annals of Applied Statistics*, **4**(1),
203-221. [doi:10.1214/09-AOAS292](https://doi.org/10.1214/09-AOAS292)

Attalides, N. (2015) Threshold-based extreme value modelling, PhD
thesis, University College London.
<https://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf>

## See also

[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md)
to set a prior density for extreme value parameters.

[`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md)
and
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md)
to perform ratio-of-uniforms sampling from an extreme value posterior
distribution.

[`kgaps_post`](https://github.com/paulnorthrop/revdbayes/reference/kgaps_post.md)
and
[`dgaps_post`](https://github.com/paulnorthrop/revdbayes/reference/dgaps_post.md)
to sample from a posterior distribution for the extremal index based on
inter-exceedance times.

The [`ru`](https://github.com/paulnorthrop/rust/reference/ru.html) and
[`ru_rcpp`](https://github.com/paulnorthrop/rust/reference/ru_rcpp.html)
functions in the
[`rust`](https://github.com/paulnorthrop/rust/reference/rust-package.html)
package for details of the arguments that can be passed to `ru` via
`rpost` and for the form of the object (of class `"evpost"`) returned
from `rpost`, which has the same structure as an object (of class
`"ru"`) returned by `ru` and `ru_rcpp`.

## Author

**Maintainer**: Paul J. Northrop <p.northrop@ucl.ac.uk> \[copyright
holder\]

Authors:

- Paul J. Northrop <p.northrop@ucl.ac.uk> \[copyright holder\]

Other contributors:

- Scott D. Grimshaw \[contributor\]
