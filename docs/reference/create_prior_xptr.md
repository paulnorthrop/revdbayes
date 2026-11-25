# Create an external pointer to a C++ prior

This function provides an example of a way in which a user can specify
their own prior density to
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md).
More specifically, a function like this (the user will need to create an
edited version tailored to their own C++ function(s)) can be used to
generate an external pointer to a compiled C++ function that evaluates
the log-prior density. Please see the vignette "Faster simulation using
revdbayes" for more information.

## Usage

``` r
create_prior_xptr(fstr)
```

## Arguments

- fstr:

  A string indicating the C++ function required.

## Value

An external pointer.

## Details

Suppose that the user's C++ functions are in a file called
"user_fns.cpp". These functions must be compiled and made available to R
before the pointer is created. This can be achieved using the function
[`sourceCpp`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) in the
**Rcpp** package or using RStudio's Source button on the editor toolbar.

For details see the examples in the documentation of the functions
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md)
and
[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md),
the vignette "Faster simulation using revdbayes" and the vignette
"Rusting Faster: Simulation using Rcpp" in the package **rust**.

## See also

[`set_prior`](https://github.com/paulnorthrop/revdbayes/reference/set_prior.md)
to specify a prior distribution using an external pointer returned by
`create_prior_xptr` and for details of in-built named prior
distributions.

The examples in the documentation of
[`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md).

## Examples

``` r
ptr_gp_flat <- create_prior_xptr("gp_flat")
prior_cfn <- set_prior(prior = ptr_gp_flat, model = "gp", min_xi = -1)

ptr_gev_flat <- create_prior_xptr("gev_flat")
prior_cfn <- set_prior(prior = ptr_gev_flat, model = "gev", min_xi = -1,
                       max_xi = Inf)

mat <- diag(c(10000, 10000, 100))
ptr_gev_norm <- create_prior_xptr("gev_norm")
pn_u <- set_prior(prior = ptr_gev_norm, model = "gev", mean = c(0,0,0),
                  icov = solve(mat))
```
