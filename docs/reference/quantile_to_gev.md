# Converts quantiles to GEV parameters

Three quantiles, that is, the value of quantile and their respective
exceedance probabilities, are provided. This function attempts to find
the location, scale and shape parameters of a GEV distribution that has
these quantiles.

## Usage

``` r
quantile_to_gev(quant, prob)
```

## Arguments

- quant:

  A numeric vector of length 3. Values of the quantiles. The values
  should *increase* with the index of the vector. If not, the values in
  `quant` will be sorted into increasing order without warning.

- prob:

  A numeric vector of length 3. Exceedance probabilities corresponding
  to the quantiles in `quant`. The values should *decrease* with the
  index of the vector. If not, the values in `prob` will be sorted into
  decreasing order without warning.

## Value

A numeric vector of length 3 containing the GEV location, scale and
shape parameters.

## Details

Suppose that \\G(x)\\ is the distribution function of a GEV(\\\mu,
\sigma, \xi\\) distribution. This function attempts to solve numerically
the set of three non-linear equations \$\$G(q_i) = 1 - p_i, i = 1, 2,
3\$\$ where \\q_i, i=1,2,3\\ are the quantiles in `quant` and \\p_i,
i=1,2,3\\ are the exceedance probabilities in `prob`. This is reduced to
a one-dimensional optimisation over the GEV shape parameter.

## See also

[`rprior_quant`](https://github.com/paulnorthrop/revdbayes/reference/rprior_quant.md)
for simulation of GEV parameters from a prior constructed on the
quantile scale.

## Examples

``` r
my_q <- c(15, 20, 22.5)
my_p <- 1-c(0.5, 0.9, 0.5^0.01)
x <- quantile_to_gev(quant = my_q, prob = my_p)
# Check
qgev(p = 1 - my_p, loc = x[1], scale = x[2], shape = x[3])
#> [1] 15.0 20.0 22.5
```
