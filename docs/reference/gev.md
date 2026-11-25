# The Generalised Extreme Value Distribution

Density function, distribution function, quantile function and random
generation for the generalised extreme value (GEV) distribution.

## Usage

``` r
dgev(x, loc = 0, scale = 1, shape = 0, log = FALSE, m = 1)

pgev(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE, m = 1)

qgev(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE, m = 1)

rgev(n, loc = 0, scale = 1, shape = 0, m = 1)
```

## Arguments

- x, q:

  Numeric vectors of quantiles.

- loc, scale, shape:

  Numeric vectors. Location, scale and shape parameters. All elements of
  `scale` must be positive.

- log, log.p:

  A logical scalar; if TRUE, probabilities p are given as log(p).

- m:

  A numeric scalar. The distribution is reparameterised by working with
  the GEV(`loc, scale, shape`) distribution function raised to the power
  `m`. See **Details**.

- lower.tail:

  A logical scalar. If TRUE (default), probabilities are \\P\[X \leq
  x\]\\, otherwise, \\P\[X \> x\]\\.

- p:

  A numeric vector of probabilities in \[0,1\].

- n:

  Numeric scalar. The number of observations to be simulated. If
  `length(n) > 1` then `length(n)` is taken to be the number required.

## Value

`dgev` gives the density function, `pgev` gives the distribution
function, `qgev` gives the quantile function, and `rgev` generates
random deviates.

The length of the result is determined by `n` for `rgev`, and is the
maximum of the lengths of the numerical arguments for the other
functions.

The numerical arguments other than `n` are recycled to the length of the
result.

## Details

The distribution function of a GEV distribution with parameters `loc` =
\\\mu\\, `scale` = \\\sigma (\> 0)\\ and `shape` = \\\xi\\ is \$\$F(x) =
\exp\\-\[1 + \xi (x - \mu) / \sigma\] ^ {-1/\xi} \\\$\$ for \\1 + \xi
(x - \mu) / \sigma \> 0\\. If \\\xi = 0\\ the distribution function is
defined as the limit as \\\xi\\ tends to zero. The support of the
distribution depends on \\\xi\\: it is \\x \leq \mu - \sigma / \xi\\ for
\\\xi \< 0\\; \\x \geq \mu - \sigma / \xi\\ for \\\xi \> 0\\; and \\x\\
is unbounded for \\\xi = 0\\. Note that if \\\xi \< -1\\ the GEV density
function becomes infinite as \\x\\ approaches \\\mu -\sigma / \xi\\ from
below.

If `lower.tail = TRUE` then if `p = 0` (`p = 1`) then the lower (upper)
limit of the distribution is returned, which is `-Inf` or `Inf` in some
cases. Similarly, but reversed, if `lower.tail = FALSE`.

See
<https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution>
for further information.

The effect of `m` is to change the location, scale and shape parameters
to \\(\mu + \sigma \log m, \sigma, \xi)\\ if \\\xi = 0\\ and \\(\mu +
\sigma (m ^ \xi - 1) / \xi, \sigma m ^ \xi, \xi)\\. For integer `m` we
can think of this as working with the maximum of `m` independent copies
of the original GEV(`loc, scale, shape`) variable.

## References

Jenkinson, A. F. (1955) The frequency distribution of the annual maximum
(or minimum) of meteorological elements. *Quart. J. R. Met. Soc.*,
**81**, 158-171.
[doi:10.1002/qj.49708134804](https://doi.org/10.1002/qj.49708134804)

Coles, S. G. (2001) *An Introduction to Statistical Modeling of Extreme
Values*, Springer-Verlag, London. Chapter 3:
[doi:10.1007/978-1-4471-3675-0_3](https://doi.org/10.1007/978-1-4471-3675-0_3)

## Examples

``` r
dgev(-1:4, 1, 0.5, 0.8)
#> [1] 0.00000000 0.00000000 0.73575888 0.17210639 0.06706381 0.03428205
dgev(1:6, 1, 0.5, -0.2, log = TRUE)
#> [1] -0.3068528 -1.4279153 -5.7449245       -Inf       -Inf       -Inf
dgev(1, shape = c(-0.2, 0.4))
#> [1] 0.2951551 0.2001168

pgev(-1:4, 1, 0.5, 0.8)
#> [1] 0.0000000 0.0000000 0.3678794 0.7386812 0.8467772 0.8948490
pgev(1:6, 1, 0.5, -0.2)
#> [1] 0.3678794 0.9251864 0.9996801 1.0000000 1.0000000 1.0000000
pgev(1, c(1, 2), c(1, 2), c(-0.2, 0.4))
#> [1] 0.3678794 0.1743086
pgev(-3, c(1, 2), c(1, 2), c(-0.2, 0.4))
#> [1] 6.218855e-09 0.000000e+00
pgev(7, 1, 1, c(-0.2, 0.4))
#> [1] 1.0000000 0.9541694

qgev((1:9)/10, 2, 0.5, 0.8)
#> [1] 1.695706 1.802111 1.913749 2.045276 2.212953 2.444700 2.800811 3.449973
#> [9] 5.157141
qgev(0.5, c(1,2), c(0.5, 1), c(-0.5, 0.5))
#> [1] 1.167445 2.402245

p <- (1:9)/10
pgev(qgev(p, 1, 2, 0.8), 1, 2, 0.8)
#> [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9

rgev(6, 1, 0.5, 0.8)
#> [1] 0.5356991 1.9642375 0.5264821 1.4918373 0.7538778 1.1022939
```
