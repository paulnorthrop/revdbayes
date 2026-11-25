# Print method for objects of class "summary.evpost"

`print` method for an object `object` of class "summary.evpost".

## Usage

``` r
# S3 method for class 'summary.evpost'
print(x, ...)
```

## Arguments

- x:

  An object of class "summary.evpost", a result of a call to
  [`summary.evpost`](https://github.com/paulnorthrop/revdbayes/reference/summary.evpost.md).

- ...:

  Additional arguments passed on to `print`.

## Value

Prints

- information about the ratio-of-uniforms bounding box, i.e.
  `object$box`

- an estimate of the probability of acceptance, i.e. `object$pa`

- a summary of the simulated values, via `summary(object$sim_vals)`

## See also

[`ru`](https://github.com/paulnorthrop/rust/reference/ru.html) or
[`ru_rcpp`](https://github.com/paulnorthrop/rust/reference/ru_rcpp.html)
for descriptions of `object$sim_vals` and `$box`.

[`plot.evpost`](https://github.com/paulnorthrop/revdbayes/reference/plot.evpost.md)
for a diagnostic plot.

## Examples

``` r
# GP posterior
u <- stats::quantile(gom, probs = 0.65)
fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
gpg <- rpost_rcpp(n = 1000, model = "gp", prior = fp, thresh = u,
                  data = gom)
summary(gpg)
#> ru bounding box:  
#>                box       vals1       vals2 conv
#> a        1.0000000  0.00000000  0.00000000    0
#> b1minus -0.1477092 -0.23028294 -0.00571716    0
#> b2minus -0.1397605  0.07992048 -0.21248500    0
#> b1plus   0.1913469  0.34032615 -0.01569524    0
#> b2plus   0.2034900  0.10510888  0.37350180    0
#> 
#> estimated probability of acceptance:  
#> [1] 0.538503
#> 
#> sample summary 
#>     sigma[u]           xi          
#>  Min.   :1.087   Min.   :-0.12897  
#>  1st Qu.:1.580   1st Qu.: 0.06139  
#>  Median :1.738   Median : 0.12524  
#>  Mean   :1.743   Mean   : 0.13165  
#>  3rd Qu.:1.900   3rd Qu.: 0.19313  
#>  Max.   :2.564   Max.   : 0.54348  
```
