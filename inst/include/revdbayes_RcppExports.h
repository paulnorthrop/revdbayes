// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_revdbayes_RCPPEXPORTS_H_GEN_
#define RCPP_revdbayes_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace revdbayes {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("revdbayes", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("revdbayes", "revdbayes_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in revdbayes");
            }
        }
    }

    inline double cpp_gp_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
        typedef SEXP(*Ptr_cpp_gp_loglik)(SEXP,SEXP);
        static Ptr_cpp_gp_loglik p_cpp_gp_loglik = NULL;
        if (p_cpp_gp_loglik == NULL) {
            validateSignature("double(*cpp_gp_loglik)(const Rcpp::NumericVector&,const Rcpp::List&)");
            p_cpp_gp_loglik = (Ptr_cpp_gp_loglik)R_GetCCallable("revdbayes", "revdbayes_cpp_gp_loglik");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpp_gp_loglik(Rcpp::wrap(x), Rcpp::wrap(ss));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline SEXP loglik_xptr(std::string fstr) {
        typedef SEXP(*Ptr_loglik_xptr)(SEXP);
        static Ptr_loglik_xptr p_loglik_xptr = NULL;
        if (p_loglik_xptr == NULL) {
            validateSignature("SEXP(*loglik_xptr)(std::string)");
            p_loglik_xptr = (Ptr_loglik_xptr)R_GetCCallable("revdbayes", "revdbayes_loglik_xptr");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_loglik_xptr(Rcpp::wrap(fstr));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline double cpp_gp_mdi(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
        typedef SEXP(*Ptr_cpp_gp_mdi)(SEXP,SEXP);
        static Ptr_cpp_gp_mdi p_cpp_gp_mdi = NULL;
        if (p_cpp_gp_mdi == NULL) {
            validateSignature("double(*cpp_gp_mdi)(const Rcpp::NumericVector&,const Rcpp::List&)");
            p_cpp_gp_mdi = (Ptr_cpp_gp_mdi)R_GetCCallable("revdbayes", "revdbayes_cpp_gp_mdi");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpp_gp_mdi(Rcpp::wrap(x), Rcpp::wrap(ppars));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double cpp_gp_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
        typedef SEXP(*Ptr_cpp_gp_flat)(SEXP,SEXP);
        static Ptr_cpp_gp_flat p_cpp_gp_flat = NULL;
        if (p_cpp_gp_flat == NULL) {
            validateSignature("double(*cpp_gp_flat)(const Rcpp::NumericVector&,const Rcpp::List&)");
            p_cpp_gp_flat = (Ptr_cpp_gp_flat)R_GetCCallable("revdbayes", "revdbayes_cpp_gp_flat");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpp_gp_flat(Rcpp::wrap(x), Rcpp::wrap(ppars));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline SEXP logprior_xptr(std::string fstr) {
        typedef SEXP(*Ptr_logprior_xptr)(SEXP);
        static Ptr_logprior_xptr p_logprior_xptr = NULL;
        if (p_logprior_xptr == NULL) {
            validateSignature("SEXP(*logprior_xptr)(std::string)");
            p_logprior_xptr = (Ptr_logprior_xptr)R_GetCCallable("revdbayes", "revdbayes_logprior_xptr");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_logprior_xptr(Rcpp::wrap(fstr));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline double cpp_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
        typedef SEXP(*Ptr_cpp_logpost)(SEXP,SEXP);
        static Ptr_cpp_logpost p_cpp_logpost = NULL;
        if (p_cpp_logpost == NULL) {
            validateSignature("double(*cpp_logpost)(const Rcpp::NumericVector&,const Rcpp::List&)");
            p_cpp_logpost = (Ptr_cpp_logpost)R_GetCCallable("revdbayes", "revdbayes_cpp_logpost");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpp_logpost(Rcpp::wrap(x), Rcpp::wrap(pars));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double cpp_logpost_phi(const Rcpp::NumericVector& phi, const Rcpp::List& pars, const SEXP& phi_to_theta_ptr) {
        typedef SEXP(*Ptr_cpp_logpost_phi)(SEXP,SEXP,SEXP);
        static Ptr_cpp_logpost_phi p_cpp_logpost_phi = NULL;
        if (p_cpp_logpost_phi == NULL) {
            validateSignature("double(*cpp_logpost_phi)(const Rcpp::NumericVector&,const Rcpp::List&,const SEXP&)");
            p_cpp_logpost_phi = (Ptr_cpp_logpost_phi)R_GetCCallable("revdbayes", "revdbayes_cpp_logpost_phi");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpp_logpost_phi(Rcpp::wrap(phi), Rcpp::wrap(pars), Rcpp::wrap(phi_to_theta_ptr));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline SEXP logpost_xptr(std::string fstr) {
        typedef SEXP(*Ptr_logpost_xptr)(SEXP);
        static Ptr_logpost_xptr p_logpost_xptr = NULL;
        if (p_logpost_xptr == NULL) {
            validateSignature("SEXP(*logpost_xptr)(std::string)");
            p_logpost_xptr = (Ptr_logpost_xptr)R_GetCCallable("revdbayes", "revdbayes_logpost_xptr");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_logpost_xptr(Rcpp::wrap(fstr));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline Rcpp::NumericVector gp_phi_to_theta(const Rcpp::NumericVector& phi, const Rcpp::List& user_args) {
        typedef SEXP(*Ptr_gp_phi_to_theta)(SEXP,SEXP);
        static Ptr_gp_phi_to_theta p_gp_phi_to_theta = NULL;
        if (p_gp_phi_to_theta == NULL) {
            validateSignature("Rcpp::NumericVector(*gp_phi_to_theta)(const Rcpp::NumericVector&,const Rcpp::List&)");
            p_gp_phi_to_theta = (Ptr_gp_phi_to_theta)R_GetCCallable("revdbayes", "revdbayes_gp_phi_to_theta");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_gp_phi_to_theta(Rcpp::wrap(phi), Rcpp::wrap(user_args));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::NumericVector >(rcpp_result_gen);
    }

    inline SEXP phi_to_theta_xptr(std::string fstr) {
        typedef SEXP(*Ptr_phi_to_theta_xptr)(SEXP);
        static Ptr_phi_to_theta_xptr p_phi_to_theta_xptr = NULL;
        if (p_phi_to_theta_xptr == NULL) {
            validateSignature("SEXP(*phi_to_theta_xptr)(std::string)");
            p_phi_to_theta_xptr = (Ptr_phi_to_theta_xptr)R_GetCCallable("revdbayes", "revdbayes_phi_to_theta_xptr");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_phi_to_theta_xptr(Rcpp::wrap(fstr));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

}

#endif // RCPP_revdbayes_RCPPEXPORTS_H_GEN_
