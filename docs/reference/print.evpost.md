# Print method for objects of class "evpost"

Print method for objects of class "evpost"

## Usage

``` r
# S3 method for class 'evpost'
print(x, ...)
```

## Arguments

- x:

  An object of class "evpost", a result of a call to
  [`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md),
  [`rpost_rcpp`](https://github.com/paulnorthrop/revdbayes/reference/rpost_rcpp.md),
  [`kgaps_post`](https://github.com/paulnorthrop/revdbayes/reference/kgaps_post.md)
  or
  [`dgaps_post`](https://github.com/paulnorthrop/revdbayes/reference/dgaps_post.md).

- ...:

  Further arguments. None are used.

## Value

The argument `x` is returned, invisibly.

## Details

`print.evpost` just prints the original function call, to avoid printing
a huge list.

## See also

[`plot.evpost`](https://github.com/paulnorthrop/revdbayes/reference/plot.evpost.md)
for a diagnostic plot.

## Examples

``` r
# Newlyn sea surges

thresh <- quantile(newlyn, probs = 0.90)
k_postsim <- kgaps_post(newlyn, thresh)
k_postsim
#> 
#> Call:
#> kgaps_post(data = newlyn, thresh = thresh)
#> 
```
