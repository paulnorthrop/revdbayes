# Simulation from a Dirichlet distribution

Simulates from a Dirichlet distribution with concentration parameter
vector \\\alpha\\ = (\\\alpha_1\\, ..., \\\alpha_K\\).

## Usage

``` r
rDir(n = 1, alpha = c(1, 1))
```

## Arguments

- n:

  A numeric scalar. The size of sample required.

- alpha:

  A numeric vector. Dirichlet concentration parameter.

## Value

An `n` by `length(alpha)` numeric matrix.

## Details

The simulation is based on the property that if \\Y_1, \ldots, Y_K\\ are
independent, \\Y_i\\ has a gamma(\\\alpha_i\\, 1) distribution and \\S =
Y_1 + \cdots + Y_k\\ then \\(Y_1, \ldots, Y_K) / S\\ has a
Dirichlet(\\\alpha_1\\, ..., \\\alpha_K\\) distribution.

See
<https://en.wikipedia.org/wiki/Dirichlet_distribution#Gamma_distribution>

## References

Kotz, S., Balakrishnan, N. and Johnson, N. L. (2000) *Continuous
Multivariate Distributions, vol. 1, Models and Applications, 2nd edn*,
ch. 49. New York: Wiley.
[doi:10.1002/0471722065](https://doi.org/10.1002/0471722065)

## See also

[`rprior_prob`](https://github.com/paulnorthrop/revdbayes/reference/rprior_prob.md)
for prior simulation of GEV parameters - prior on probability scale.

## Examples

``` r
rDir(n = 10, alpha = 1:4)
#>              [,1]       [,2]       [,3]      [,4]
#>  [1,] 0.115370906 0.01462583 0.23753035 0.6324729
#>  [2,] 0.049415940 0.49527651 0.09582583 0.3594817
#>  [3,] 0.028259268 0.11340945 0.24134632 0.6169850
#>  [4,] 0.123849100 0.23267500 0.32954730 0.3139286
#>  [5,] 0.063862382 0.12346098 0.17856201 0.6341146
#>  [6,] 0.079358460 0.14770560 0.25000066 0.5229353
#>  [7,] 0.001612992 0.24626422 0.07988016 0.6722426
#>  [8,] 0.080030390 0.02010660 0.45993851 0.4399245
#>  [9,] 0.023265218 0.29859926 0.42866548 0.2494700
#> [10,] 0.037039164 0.02016963 0.23130954 0.7114817
```
