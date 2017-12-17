#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    _revdbayes_any_nonpos
    _revdbayes_cpp_gp_loglik
    _revdbayes_cpp_gev_loglik
    _revdbayes_cpp_os_loglik
    _revdbayes_cpp_pp_loglik
    _revdbayes_cpp_gp_norm
    _revdbayes_cpp_gp_mdi
    _revdbayes_cpp_gp_flat
    _revdbayes_cpp_gp_flatflat
    _revdbayes_cpp_gp_jeffreys
    _revdbayes_cpp_gp_beta
    _revdbayes_cpp_gev_norm
    _revdbayes_cpp_gev_loglognorm
    _revdbayes_cpp_gev_mdi
    _revdbayes_cpp_gev_flat
    _revdbayes_cpp_gev_flatflat
    _revdbayes_cpp_gev_beta
    _revdbayes_lgdgev_cpp
    _revdbayes_pgev_cpp
    _revdbayes_qgev_cpp
    _revdbayes_cpp_gev_prob
    _revdbayes_cpp_gev_quant
    _revdbayes_gp_user_logpost
    _revdbayes_gev_user_logpost
    _revdbayes_os_user_logpost
    _revdbayes_pp_user_logpost
    _revdbayes_gp_mdi_logpost
    _revdbayes_gp_norm_logpost
    _revdbayes_gp_flat_logpost
    _revdbayes_gp_flatflat_logpost
    _revdbayes_gp_jeffreys_logpost
    _revdbayes_gp_beta_logpost
    _revdbayes_gev_mdi_logpost
    _revdbayes_gev_norm_logpost
    _revdbayes_gev_loglognorm_logpost
    _revdbayes_gev_flat_logpost
    _revdbayes_gev_flatflat_logpost
    _revdbayes_gev_beta_logpost
    _revdbayes_gev_prob_logpost
    _revdbayes_gev_quant_logpost
    _revdbayes_pp_mdi_logpost
    _revdbayes_pp_norm_logpost
    _revdbayes_pp_loglognorm_logpost
    _revdbayes_pp_flat_logpost
    _revdbayes_pp_flatflat_logpost
    _revdbayes_pp_beta_logpost
    _revdbayes_pp_prob_logpost
    _revdbayes_pp_quant_logpost
    _revdbayes_os_mdi_logpost
    _revdbayes_os_norm_logpost
    _revdbayes_os_loglognorm_logpost
    _revdbayes_os_flat_logpost
    _revdbayes_os_flatflat_logpost
    _revdbayes_os_beta_logpost
    _revdbayes_os_prob_logpost
    _revdbayes_os_quant_logpost
    _revdbayes_gp_logpost_xptr
    _revdbayes_gev_logpost_xptr
    _revdbayes_os_logpost_xptr
    _revdbayes_pp_logpost_xptr
    _revdbayes_gp_phi_to_theta
    _revdbayes_gev_phi_to_theta
    _revdbayes_pp_phi_to_theta
    _revdbayes_kgaps_phi_to_theta
    _revdbayes_phi_to_theta_xptr
    _revdbayes_gp_mdi_logpost_phi
    _revdbayes_gp_norm_logpost_phi
    _revdbayes_gp_flat_logpost_phi
    _revdbayes_gp_flatflat_logpost_phi
    _revdbayes_gp_jeffreys_logpost_phi
    _revdbayes_gp_beta_logpost_phi
    _revdbayes_gp_user_logpost_phi
    _revdbayes_gev_mdi_logpost_phi
    _revdbayes_gev_norm_logpost_phi
    _revdbayes_gev_loglognorm_logpost_phi
    _revdbayes_gev_flat_logpost_phi
    _revdbayes_gev_flatflat_logpost_phi
    _revdbayes_gev_beta_logpost_phi
    _revdbayes_gev_prob_logpost_phi
    _revdbayes_gev_quant_logpost_phi
    _revdbayes_gev_user_logpost_phi
    _revdbayes_os_mdi_logpost_phi
    _revdbayes_os_norm_logpost_phi
    _revdbayes_os_loglognorm_logpost_phi
    _revdbayes_os_flat_logpost_phi
    _revdbayes_os_flatflat_logpost_phi
    _revdbayes_os_beta_logpost_phi
    _revdbayes_os_prob_logpost_phi
    _revdbayes_os_quant_logpost_phi
    _revdbayes_os_user_logpost_phi
    _revdbayes_pp_mdi_logpost_phi
    _revdbayes_pp_norm_logpost_phi
    _revdbayes_pp_loglognorm_logpost_phi
    _revdbayes_pp_flat_logpost_phi
    _revdbayes_pp_flatflat_logpost_phi
    _revdbayes_pp_beta_logpost_phi
    _revdbayes_pp_prob_logpost_phi
    _revdbayes_pp_quant_logpost_phi
    _revdbayes_pp_user_logpost_phi
    _revdbayes_gp_logpost_phi_xptr
    _revdbayes_gev_logpost_phi_xptr
    _revdbayes_pp_logpost_phi_xptr
    _revdbayes_os_logpost_phi_xptr
    _revdbayes_kgaps_logpost
    _revdbayes_kgaps_logpost_xptr
    _revdbayes_kgaps_log_j
    _revdbayes_log_j_xptr
    _revdbayes_user_gp_flat
    _revdbayes_user_gev_norm
    _revdbayes_user_gev_flat
    _revdbayes_create_prior_xptr

  Most likely possible values need to be added below.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _revdbayes_RcppExport_registerCCallable();

static const R_CallMethodDef CallEntries[] = {
    {"_revdbayes_RcppExport_registerCCallable", (DL_FUNC) &_revdbayes_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

void R_init_revdbayes(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
