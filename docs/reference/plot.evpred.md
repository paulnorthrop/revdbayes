# Plot diagnostics for an evpred object

`plot` method for class "evpred". Plots summarising the predictive
distribution of the largest value to be observed in N years are
produced. The plot produced depends on `x$type`. If `x$type = "d", "p"`
or `"q"` then [`matplot`](https://rdrr.io/r/graphics/matplot.html) is
used to produce a line plot of the predictive density, distribution or
quantile function, respectively, with a line for each value of N in
`x$n_years`. If `x$type = "r"` then estimates of the predictive density
(from [`density`](https://rdrr.io/r/stats/density.html)) are plotted
with a line for each N. If `x$type = "i"` then lines representing
estimated predictive intervals are plotted, with the level of the
interval indicated next to the line.

## Usage

``` r
# S3 method for class 'evpred'
plot(
  x,
  ...,
  leg_pos = NULL,
  leg_text = NULL,
  which_int = c("long", "short", "both")
)
```

## Arguments

- x:

  An object of class "evpost", a result of a call to
  [`rpost`](https://github.com/paulnorthrop/revdbayes/reference/rpost.md).

- ...:

  Additional arguments passed on to
  [`matplot`](https://rdrr.io/r/graphics/matplot.html).

- leg_pos:

  A character scalar. Keyword for the position of legend. See
  [`legend`](https://rdrr.io/r/graphics/legend.html).

- leg_text:

  A character or expression vector. Text for legend. See
  [`legend`](https://rdrr.io/r/graphics/legend.html).

- which_int:

  A character scalar. If `x$type = "i"` which intervals should be
  plotted? `"long"` for equi-tailed intervals, `"short"` for the
  shortest possible intervals, `"both"` for both.

## Value

Nothing is returned.

## See also

[`predict.evpost`](https://github.com/paulnorthrop/revdbayes/reference/predict.evpost.md)
for the S3 `predict` method for objects of class `evpost`.

## Examples

``` r
data(portpirie)
mat <- diag(c(10000, 10000, 100))
pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
gevp  <- rpost(n = 1000, model = "gev", prior = pn, data = portpirie)

# Predictive density function
d_gevp <- predict(gevp, type = "d", n_years = c(100, 1000))
plot(d_gevp)


# Predictive distribution function
p_gevp <- predict(gevp, type = "p", n_years = c(100, 1000))
plot(p_gevp)


# Predictive quantiles
q_gevp <- predict(gevp, type = "q", n_years = c(100, 1000))
plot(q_gevp)


# Predictive intervals
i_gevp <- predict(gevp, type = "i", n_years = c(100, 1000), hpd = TRUE)
plot(i_gevp, which_int = "both")

#> Warning: argument 1 does not name a graphical parameter

# Sample from predictive distribution
r_gevp <- predict(gevp, type = "r", n_years = c(100, 1000))
plot(r_gevp)

plot(r_gevp, xlim = c(4, 10))

```
