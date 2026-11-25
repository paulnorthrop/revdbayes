# The Generalised Pareto Distribution

Density function, distribution function, quantile function and random
generation for the generalised Pareto (GP) distribution.

## Usage

``` r
dgp(x, loc = 0, scale = 1, shape = 0, log = FALSE)

pgp(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)

qgp(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)

rgp(n, loc = 0, scale = 1, shape = 0)
```

## Arguments

- x, q:

  Numeric vectors of quantiles. All elements of `x` and `q` must be
  non-negative.

- loc, scale, shape:

  Numeric vectors. Location, scale and shape parameters. All elements of
  `scale` must be positive.

- log, log.p:

  A logical scalar; if TRUE, probabilities p are given as log(p).

- lower.tail:

  A logical scalar. If TRUE (default), probabilities are \\P\[X \leq
  x\]\\, otherwise, \\P\[X \> x\]\\.

- p:

  A numeric vector of probabilities in \[0,1\].

- n:

  Numeric scalar. The number of observations to be simulated. If
  `length(n) > 1` then `length(n)` is taken to be the number required.

## Value

`dgp` gives the density function, `pgp` gives the distribution function,
`qgp` gives the quantile function, and `rgp` generates random deviates.

## Details

The distribution function of a GP distribution with parameters
`location` = \\\mu\\, `scale` = \\\sigma (\> 0)\\ and `shape` = \\\xi\\
is \$\$F(x) = 1 - \[1 + \xi (x - \mu) / \sigma\] ^ {-1/\xi}\$\$ for
\\1 + \xi (x - \mu) / \sigma \> 0\\. If \\\xi = 0\\ the distribution
function is defined as the limit as \\\xi\\ tends to zero. The support
of the distribution depends on \\\xi\\: it is \\x \geq \mu\\ for \\\xi
\geq 0\\; and \\\mu \leq x \leq \mu - \sigma / \xi\\ for \\\xi \< 0\\.
Note that if \\\xi \< -1\\ the GP density function becomes infinite as
\\x\\ approaches \\\mu - \sigma/\xi\\.

If `lower.tail = TRUE` then if `p = 0` (`p = 1`) then the lower (upper)
limit of the distribution is returned. The upper limit is `Inf` if
`shape` is non-negative. Similarly, but reversed, if
`lower.tail = FALSE`.

See <https://en.wikipedia.org/wiki/Generalized_Pareto_distribution> for
further information.

## References

Pickands, J. (1975) Statistical inference using extreme order
statistics. *Annals of Statistics*, **3**, 119-131.
[doi:10.1214/aos/1176343003](https://doi.org/10.1214/aos/1176343003)

Coles, S. G. (2001) *An Introduction to Statistical Modeling of Extreme
Values*, Springer-Verlag, London. Chapter 4:
[doi:10.1007/978-1-4471-3675-0_4](https://doi.org/10.1007/978-1-4471-3675-0_4)

## Examples

``` r
dgp(0:4, scale = 0.5, shape = 0.8)
#> [1] 2.00000000 0.23299144 0.07919889 0.03831043 0.02214413
dgp(1:6, scale = 0.5, shape = -0.2, log = TRUE)
#> [1] -1.350155 -5.744604      -Inf      -Inf      -Inf      -Inf
dgp(1, scale = 1, shape = c(-0.2, 0.4))
#> [1] 0.4096000 0.3080008

pgp(0:4, scale = 0.5, shape = 0.8)
#> [1] 0.0000000 0.6971111 0.8336823 0.8888998 0.9180667
pgp(1:6, scale = 0.5, shape = -0.2)
#> [1] 0.92224 0.99968 1.00000 1.00000 1.00000 1.00000
pgp(1, scale = c(1, 2), shape = c(-0.2, 0.4))
#> [1] 0.6723200 0.3660619
pgp(7, scale = 1, shape = c(-0.2, 0.4))
#> [1] 1.0000000 0.9644744

qgp((0:9)/10, scale = 0.5, shape = 0.8)
#>  [1] 0.00000000 0.05496414 0.12215039 0.20638385 0.31550047 0.46318820
#>  [7] 0.67586439 1.01250643 1.63993645 3.31848340
qgp(0.5, scale = c(0.5, 1), shape = c(-0.5, 0.5))
#> [1] 0.2928932 0.8284271

p <- (1:9)/10
pgp(qgp(p, scale = 2, shape = 0.8), scale = 2, shape = 0.8)
#> [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9

rgp(6, scale = 0.5, shape = 0.8)
#> [1] 1.2688250 0.4340321 2.5833754 0.3176955 0.1650756 0.4685481
```
