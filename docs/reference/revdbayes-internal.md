# Internal revdbayes functions

Internal revdbayes functions

## Usage

``` r
process_data(model, data, thresh, noy, use_noy, ros, weights = NULL)

create_ru_list(model, trans, rotate, min_xi, max_xi)

set_which_lam(model)

set_range_phi(model, phi_mid, se_phi, mult)

box_cox(x, lambda = 1, gm = 1, lambda_tol = 1e-06, poly_order = 3)

box_cox_vec(x, lambda = 1, lambda_tol = 1e-06)

box_cox_deriv(x, lambda = 1, lambda_tol = 1e-06, poly_order = 3)

check_sample_size(prior_name, n_check)

check_sample_size_message(prior_name, n_check)

logNegNA(x)

gp_mle(gp_data)

fallback_gp_mle(init, ...)
```

## Details

These functions are not intended to be called by the user.
