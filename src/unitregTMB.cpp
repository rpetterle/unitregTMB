// ----------------------------------------------------------------------------
// unitregTMB.cpp
// ----------------------------------------------------------------------------

#define EIGEN_DONT_PARALLELIZE
#define TMB_LIB_INIT tmb_custom_init

#include <TMB.hpp>
#include <cmath>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "Distributions.h"


enum valid_mu_link {
  logit_link   = 0,
  probit_link  = 1,
  cloglog_link = 2,
  cauchit_link = 3
};


enum valid_family {
  // Mean
  fam_beta_mean          = 0,
  fam_simplex_mean       = 1,
  fam_vasicek_mean       = 2,
  fam_ugamma_mean        = 3,
  fam_bessel_mean        = 4,

  // Mode
  fam_beta_mode          = 5,
  fam_kum_mode           = 6,
  fam_ugamma_mode        = 7,
  fam_ugompertz_mode     = 8,

  // Quantile
  fam_kum_quantile       = 9,
  fam_vasicek_quantile   = 10,
  fam_uweibull_quantile  = 11,
  fam_ugompertz_quantile = 12,
  fam_johnsonsb_quantile = 13,
  fam_ashw_quantile      = 14,
  fam_ubs_quantile       = 15
};


template<class Type>
Type inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
    case logit_link:
      ans = invlogit(eta);
      break;
    case probit_link:
      ans = pnorm(eta);
      break;
    case cloglog_link:
      ans = Type(1.0) - exp(-exp(eta));
      break;
    case cauchit_link:
      ans = Type(0.5) + atan(eta) / Type(M_PI);
      break;
    default:
      error("Link not implemented or invalid for mu!");
  }
  return ans;
}


template<class Type>
Type get_log_density(int family, Type y, Type mu, Type phi, Type tau) {
  switch (family) {
    // Mean
    case fam_beta_mean:          return dbeta_mean(y, mu, phi, 1);
    case fam_simplex_mean:       return dsimplex(y, mu, phi, 1);
    case fam_vasicek_mean:       return dvasicek_mean(y, mu, phi, 1);
    case fam_ugamma_mean:        return dunitgamma_mean(y, mu, phi, 1);
    case fam_bessel_mean:        return dbessel_mean(y, mu, phi, 1);

    // Mode
    case fam_beta_mode:          return dbeta_mode(y, mu, phi, 1);
    case fam_kum_mode:           return dkum_mode(y, mu, phi, 1);
    case fam_ugamma_mode:        return dugamma_mode(y, mu, phi, 1);
    case fam_ugompertz_mode:     return dunitgompertz_mode(y, mu, phi, 1);

    // Quantile
    case fam_kum_quantile:       return dkum_quantile(y, mu, phi, tau, 1);
    case fam_vasicek_quantile:   return dvasicek_quantile(y, mu, phi, tau, 1);
    case fam_uweibull_quantile:  return dunitweibull(y, mu, phi, tau, 1);
    case fam_ugompertz_quantile: return dunitgompertz_quantile(y, mu, phi, tau, 1);
    case fam_johnsonsb_quantile: return djohnsonsb(y, mu, phi, tau, 1);
    case fam_ashw_quantile:      return dashw(y, mu, phi, tau, 1);
    case fam_ubs_quantile:       return dubs(y, mu, phi, tau, 1);

    default:
      error("Unsupported family code.");
      return Type(0.0);
  }
}


template<class Type>
Type inverse_phi_linkfun(Type eta, int family) {
  switch (family) {
    case fam_vasicek_mean:
    case fam_vasicek_quantile:
      return invlogit(eta);
    default:
      return exp(eta);
  }
}

enum re_covariance_type {
    re_cov_unstructured = 0
};

inline SEXP unitreg_get_list_element(SEXP x, const char* name) {
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (Rf_isNull(names)) return R_NilValue;

    for (int i = 0; i < Rf_length(x); ++i) {
        if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            return VECTOR_ELT(x, i);
        }
    }
    return R_NilValue;
}

inline int unitreg_get_int_element(SEXP x, const char* name) {
    SEXP ans = unitreg_get_list_element(x, name);
    if (Rf_isNull(ans)) {
        error("Invalid random-effects term metadata.");
    }
    return Rf_asInteger(ans);
}

template <class Type>
struct re_term_info {
    int cov_type;
    int n_coef;
    int n_levels;
    int n_cov_par;
    int full_corr;
};

template <class Type>
struct re_terms_t : vector<re_term_info<Type> > {
    re_terms_t(SEXP x) {
        (*this).resize(Rf_length(x));

        for (int i = 0; i < Rf_length(x); ++i) {
            SEXP term = VECTOR_ELT(x, i);

            (*this)(i).cov_type =
                unitreg_get_int_element(term, "cov_type");
            (*this)(i).n_coef =
                unitreg_get_int_element(term, "n_coef");
            (*this)(i).n_levels =
                unitreg_get_int_element(term, "n_levels");
            (*this)(i).n_cov_par =
                unitreg_get_int_element(term, "n_cov_par");
            (*this)(i).full_corr =
                unitreg_get_int_element(term, "full_corr");
        }
    }
};

template <class Type>
inline void get_re_dimensions(
    const vector<re_term_info<Type> >& terms,
    int& total_n_comp,
    int& total_n_cor
) {
    total_n_comp = 0;
    total_n_cor = 0;

    for (int t = 0; t < terms.size(); ++t) {
        const int nc = terms(t).n_coef;
        total_n_comp += nc;
        total_n_cor += nc * (nc - 1) / 2;
    }
}


template <class Type>
Type re_term_nll(
    array<Type>& U,
    const vector<Type>& theta,
    const re_term_info<Type>& term,
    vector<Type>& re_sd,
    vector<Type>& re_corr,
    int& idx_std,
    int& idx_cor
) {
    using namespace density;

    if (term.cov_type != re_cov_unstructured) {
        error("Random-effects covariance structure not implemented.");
    }

    const int n_comp = term.n_coef;
    const int expected_theta = n_comp * (n_comp + 1) / 2;

    if (term.n_cov_par != expected_theta ||
        theta.size() != expected_theta) {
        error("Inconsistent random-effects covariance parameter count.");
    }

    matrix<Type> L(n_comp, n_comp);
    L.setZero();

    int k_param = 0;
    for (int j = 0; j < n_comp; ++j) {
        for (int i = j; i < n_comp; ++i) {
            if (i == j) {
                L(i, j) = exp(theta(k_param));
            } else {
                L(i, j) = theta(k_param);
            }
            ++k_param;
        }
    }

    matrix<Type> Sigma = L * L.transpose();

    for (int j = 0; j < n_comp; ++j) {
        re_sd(idx_std++) = sqrt(Sigma(j, j));

        for (int i = j + 1; i < n_comp; ++i) {
            re_corr(idx_cor++) =
                Sigma(i, j) / sqrt(Sigma(i, i) * Sigma(j, j));
        }
    }

    MVNORM_t<Type> neg_log_dmvnorm(Sigma);
    Type ans = Type(0.0);

    for (int k = 0; k < term.n_levels; ++k) {
        ans += neg_log_dmvnorm(U.col(k));
    }

    return ans;
}


template <class Type>
Type re_component_nll(
    vector<Type>& u,
    vector<Type>& theta,
    vector<re_term_info<Type> >& terms,
    vector<Type>& re_sd,
    vector<Type>& re_corr
) {
    Type ans = Type(0.0);

    int u_pointer = 0;
    int theta_pointer = 0;
    int idx_std = 0;
    int idx_cor = 0;

    for (int t = 0; t < terms.size(); ++t) {
        const int nr = terms(t).n_coef * terms(t).n_levels;
        const int np = terms(t).n_cov_par;

        vector<int> dim(2);
        dim << terms(t).n_coef, terms(t).n_levels;

        array<Type> U(&u(u_pointer), dim);
        vector<Type> theta_term = theta.segment(theta_pointer, np);

        ans += re_term_nll(
            U,
            theta_term,
            terms(t),
            re_sd,
            re_corr,
            idx_std,
            idx_cor
        );

        u_pointer += nr;
        theta_pointer += np;
    }

    if (u_pointer != u.size() || theta_pointer != theta.size()) {
        error("Random-effects metadata do not match parameter-vector sizes.");
    }

    return ans;
}


template<class Type>
Type process_spatial_effects(
    const vector<Type>& u_spatial,
    const vector<Type>& spatial_par,
    const Eigen::SparseMatrix<Type>& W_spatial,
    const Eigen::SparseMatrix<Type>& D_spatial,
    Type& spatial_scale,
    Type& spatial_rho
) {
    using namespace density;
    spatial_scale = exp(spatial_par(0));
    Type rho_upper = Type(1.0) - Type(1e-8);
    spatial_rho = rho_upper * invlogit(spatial_par(1));

    Eigen::SparseMatrix<Type> Q_spatial = D_spatial - spatial_rho * W_spatial;

    return SCALE(GMRF(Q_spatial), spatial_scale)(u_spatial);
}

// --------------------------------------------------------------------------
// OBJECTIVE FUNCTION
// --------------------------------------------------------------------------

template<class Type>
Type objective_function<Type>::operator()() {

// --------------------------------------------------------------------------
// DATA
// --------------------------------------------------------------------------

  DATA_VECTOR(Y);
  DATA_IVECTOR(y_class);

  DATA_MATRIX(X_mu);
  DATA_MATRIX(X_phi);
  DATA_MATRIX(X_p0);
  DATA_MATRIX(X_p1);

  DATA_SPARSE_MATRIX(Z_mu);
  DATA_SPARSE_MATRIX(Z_phi);
  DATA_SPARSE_MATRIX(Z_p0);
  DATA_SPARSE_MATRIX(Z_p1);

  DATA_VECTOR(weights);
  DATA_VECTOR(offset_mu);

  DATA_INTEGER(family);
  DATA_INTEGER(link);

  DATA_INTEGER(has_p0_inflation);
  DATA_INTEGER(has_p1_inflation);

  DATA_SCALAR(tau);

  DATA_STRUCT(re_terms_mu, re_terms_t);
  DATA_STRUCT(re_terms_phi, re_terms_t);
  DATA_STRUCT(re_terms_p0, re_terms_t);
  DATA_STRUCT(re_terms_p1, re_terms_t);

  const bool has_random_effects_mu  = (re_terms_mu.size()  > 0);
  const bool has_random_effects_phi = (re_terms_phi.size() > 0);
  const bool has_random_effects_p0  = (re_terms_p0.size()  > 0);
  const bool has_random_effects_p1  = (re_terms_p1.size()  > 0);

  DATA_INTEGER(has_spatial_effects_mu);
  DATA_SPARSE_MATRIX(Z_spatial_mu);
  DATA_SPARSE_MATRIX(W_spatial_mu);
  DATA_VECTOR(spatial_num_neighbors_mu);  
  DATA_SPARSE_MATRIX(D_spatial_mu);

  DATA_INTEGER(has_spatial_effects_phi);
  DATA_SPARSE_MATRIX(Z_spatial_phi);
  DATA_SPARSE_MATRIX(W_spatial_phi);
  DATA_VECTOR(spatial_num_neighbors_phi);  
  DATA_SPARSE_MATRIX(D_spatial_phi);

  DATA_INTEGER(has_spatial_effects_p0);
  DATA_SPARSE_MATRIX(Z_spatial_p0);
  DATA_SPARSE_MATRIX(W_spatial_p0);
  DATA_VECTOR(spatial_num_neighbors_p0);  
  DATA_SPARSE_MATRIX(D_spatial_p0);

  DATA_INTEGER(has_spatial_effects_p1);
  DATA_SPARSE_MATRIX(Z_spatial_p1);
  DATA_SPARSE_MATRIX(W_spatial_p1);
  DATA_VECTOR(spatial_num_neighbors_p1);  
  DATA_SPARSE_MATRIX(D_spatial_p1);

  // --------------------------------------------------------------------------
  // PARAMETERS
  // --------------------------------------------------------------------------

  PARAMETER_VECTOR(beta_mu);
  PARAMETER_VECTOR(beta_phi);
  PARAMETER_VECTOR(beta_p0);
  PARAMETER_VECTOR(beta_p1);

  PARAMETER_VECTOR(re_coef_mu);
  PARAMETER_VECTOR(re_chol_par_mu);

  PARAMETER_VECTOR(re_coef_phi);
  PARAMETER_VECTOR(re_chol_par_phi);

  PARAMETER_VECTOR(re_coef_p0);
  PARAMETER_VECTOR(re_chol_par_p0);

  PARAMETER_VECTOR(re_coef_p1);
  PARAMETER_VECTOR(re_chol_par_p1);

  PARAMETER_VECTOR(u_spatial_mu);
  PARAMETER_VECTOR(spatial_par_mu);

  PARAMETER_VECTOR(u_spatial_phi);
  PARAMETER_VECTOR(spatial_par_phi);

  PARAMETER_VECTOR(u_spatial_p0);
  PARAMETER_VECTOR(spatial_par_p0);

  PARAMETER_VECTOR(u_spatial_p1);
  PARAMETER_VECTOR(spatial_par_p1);

 // --------------------------------------------------------------------------
 // INITIALIZATION
 // --------------------------------------------------------------------------

  Type nll = Type(0.0);
  int n_obs = Y.size();

  vector<Type> conditional_log_lik_unweighted(n_obs);
  vector<Type> conditional_log_lik_weighted(n_obs);

  conditional_log_lik_unweighted.setZero();
  conditional_log_lik_weighted.setZero();


 // --------------------------------------------------------------------------
 // LINEAR PREDICTORS
 // --------------------------------------------------------------------------

  vector<Type> eta_mu = X_mu * beta_mu + offset_mu;
  if (has_random_effects_mu == 1) {
    eta_mu += Z_mu * re_coef_mu;
  }
  if (has_spatial_effects_mu == 1) {
    eta_mu += Z_spatial_mu * u_spatial_mu;
  }

  vector<Type> eta_phi = X_phi * beta_phi;
  if (has_random_effects_phi == 1) {
    eta_phi += Z_phi * re_coef_phi;
  }
  if (has_spatial_effects_phi == 1) {
    eta_phi += Z_spatial_phi * u_spatial_phi;
  }

  vector<Type> eta_p0(n_obs);
  eta_p0.setZero();
  if (has_p0_inflation == 1) {
    eta_p0 = X_p0 * beta_p0;
    if (has_random_effects_p0 == 1) {
      eta_p0 += Z_p0 * re_coef_p0;
    }
    if (has_spatial_effects_p0 == 1) {
      eta_p0 += Z_spatial_p0 * u_spatial_p0;
    }
  }

  vector<Type> eta_p1(n_obs);
  eta_p1.setZero();
  if (has_p1_inflation == 1) {
    eta_p1 = X_p1 * beta_p1;
    if (has_random_effects_p1 == 1) {
      eta_p1 += Z_p1 * re_coef_p1;
    }
    if (has_spatial_effects_p1 == 1) {
      eta_p1 += Z_spatial_p1 * u_spatial_p1;
    }
  }


 // --------------------------------------------------------------------------
 // RANDOM EFFECTS
 // --------------------------------------------------------------------------

  int total_n_comp_mu, total_n_cor_mu;
  get_re_dimensions(re_terms_mu, total_n_comp_mu, total_n_cor_mu);
  vector<Type> re_sd_mu(total_n_comp_mu);
  vector<Type> re_corr_mu(total_n_cor_mu > 0 ? total_n_cor_mu : 1);
  if (has_random_effects_mu) {
    nll += re_component_nll(re_coef_mu, re_chol_par_mu, re_terms_mu,
                            re_sd_mu, re_corr_mu);
  }

  int total_n_comp_phi, total_n_cor_phi;
  get_re_dimensions(re_terms_phi, total_n_comp_phi, total_n_cor_phi);
  vector<Type> re_sd_phi(total_n_comp_phi);
  vector<Type> re_corr_phi(total_n_cor_phi > 0 ? total_n_cor_phi : 1);
  if (has_random_effects_phi) {
    nll += re_component_nll(re_coef_phi, re_chol_par_phi, re_terms_phi,
                            re_sd_phi, re_corr_phi);
  }

  int total_n_comp_p0, total_n_cor_p0;
  get_re_dimensions(re_terms_p0, total_n_comp_p0, total_n_cor_p0);
  vector<Type> re_sd_p0(total_n_comp_p0);
  vector<Type> re_corr_p0(total_n_cor_p0 > 0 ? total_n_cor_p0 : 1);
  if (has_random_effects_p0) {
    nll += re_component_nll(re_coef_p0, re_chol_par_p0, re_terms_p0,
                            re_sd_p0, re_corr_p0);
  }

  int total_n_comp_p1, total_n_cor_p1;
  get_re_dimensions(re_terms_p1, total_n_comp_p1, total_n_cor_p1);
  vector<Type> re_sd_p1(total_n_comp_p1);
  vector<Type> re_corr_p1(total_n_cor_p1 > 0 ? total_n_cor_p1 : 1);
  if (has_random_effects_p1) {
    nll += re_component_nll(re_coef_p1, re_chol_par_p1, re_terms_p1,
                            re_sd_p1, re_corr_p1);
  }

 // --------------------------------------------------------------------------
 // SPATIAL RANDOM EFFECTS 
 // --------------------------------------------------------------------------

  Type spatial_scale_mu = 0, spatial_rho_mu = 0;
  if (has_spatial_effects_mu == 1) {
    nll += process_spatial_effects(
      u_spatial_mu, spatial_par_mu, W_spatial_mu, D_spatial_mu,
      spatial_scale_mu, spatial_rho_mu);
  }

  Type spatial_scale_phi = 0, spatial_rho_phi = 0;
  if (has_spatial_effects_phi == 1) {
    nll += process_spatial_effects(
      u_spatial_phi, spatial_par_phi, W_spatial_phi, D_spatial_phi,
      spatial_scale_phi, spatial_rho_phi);
  }

  Type spatial_scale_p0 = 0, spatial_rho_p0 = 0;
  if (has_spatial_effects_p0 == 1) {
    nll += process_spatial_effects(
      u_spatial_p0, spatial_par_p0, W_spatial_p0, D_spatial_p0,
      spatial_scale_p0, spatial_rho_p0);
  }

  Type spatial_scale_p1 = 0, spatial_rho_p1 = 0;
  if (has_spatial_effects_p1 == 1) {
    nll += process_spatial_effects(
      u_spatial_p1, spatial_par_p1, W_spatial_p1, D_spatial_p1,
      spatial_scale_p1, spatial_rho_p1);
  }

 // --------------------------------------------------------------------------
 // CONDITIONAL LIKELIHOOD
 // --------------------------------------------------------------------------

  for (int i = 0; i < n_obs; ++i) {
    Type mu_i = inverse_linkfun(eta_mu(i), link);
    Type phi_i = inverse_phi_linkfun(eta_phi(i), family);

    Type log_p0 = Type(0.0);
    Type log_p1 = Type(0.0);
    Type log_pc = Type(0.0);
    Type li = Type(0.0);

    if (has_p0_inflation == 1 && has_p1_inflation == 1) {
      Type log_boundary_sum = logspace_add(eta_p0(i), eta_p1(i));
      Type log_denom = logspace_add(Type(0.0), log_boundary_sum);
      log_p0 = eta_p0(i) - log_denom;
      log_p1 = eta_p1(i) - log_denom;
      log_pc = -log_denom;
    } else if (has_p0_inflation == 1) {
      Type log_denom = logspace_add(Type(0.0), eta_p0(i));
      log_p0 = eta_p0(i) - log_denom;
      log_pc = -log_denom;
    } else if (has_p1_inflation == 1) {
      Type log_denom = logspace_add(Type(0.0), eta_p1(i));
      log_p1 = eta_p1(i) - log_denom;
      log_pc = -log_denom;
    }

    if (y_class(i) == 0) {
      li = log_p0;
    } else if (y_class(i) == 1) {
      li = log_p1;
    } else if (y_class(i) == 2) {
      Type dens = get_log_density(family, Y(i), mu_i, phi_i, tau);
      li = log_pc + dens;
    } else {
      error("Invalid response classification.");
    }

    conditional_log_lik_unweighted(i) = li;
    li *= weights(i);
    conditional_log_lik_weighted(i) = li;
    nll -= li;
  }

 // --------------------------------------------------------------------------
 // REPORT
 // --------------------------------------------------------------------------

  REPORT(conditional_log_lik_unweighted);
  REPORT(conditional_log_lik_weighted);

  if (has_random_effects_mu == 1) {
    REPORT(re_sd_mu);
    if (total_n_cor_mu > 0) { REPORT(re_corr_mu); }
    ADREPORT(re_chol_par_mu);
    ADREPORT(re_sd_mu);
    if (total_n_cor_mu > 0) { ADREPORT(re_corr_mu); }
  }

  if (has_random_effects_phi == 1) {
    REPORT(re_sd_phi);
    if (total_n_cor_phi > 0) { REPORT(re_corr_phi); }
    ADREPORT(re_chol_par_phi);
    ADREPORT(re_sd_phi);
    if (total_n_cor_phi > 0) { ADREPORT(re_corr_phi); }
  }

  if (has_random_effects_p0 == 1) {
    REPORT(re_sd_p0);
    if (total_n_cor_p0 > 0) { REPORT(re_corr_p0); }
    ADREPORT(re_chol_par_p0);
    ADREPORT(re_sd_p0);
    if (total_n_cor_p0 > 0) { ADREPORT(re_corr_p0); }
  }

  if (has_random_effects_p1 == 1) {
    REPORT(re_sd_p1);
    if (total_n_cor_p1 > 0) { REPORT(re_corr_p1); }
    ADREPORT(re_chol_par_p1);
    ADREPORT(re_sd_p1);
    if (total_n_cor_p1 > 0) { ADREPORT(re_corr_p1); }
  }

  if (has_spatial_effects_mu == 1) {
    REPORT(spatial_scale_mu);
    REPORT(spatial_rho_mu);
    REPORT(u_spatial_mu);
    ADREPORT(spatial_scale_mu);
    ADREPORT(spatial_rho_mu);
  }

  if (has_spatial_effects_phi == 1) {
    REPORT(spatial_scale_phi);
    REPORT(spatial_rho_phi);
    REPORT(u_spatial_phi);
    ADREPORT(spatial_scale_phi);
    ADREPORT(spatial_rho_phi);
  }

  if (has_spatial_effects_p0 == 1) {
    REPORT(spatial_scale_p0);
    REPORT(spatial_rho_p0);
    REPORT(u_spatial_p0);
    ADREPORT(spatial_scale_p0);
    ADREPORT(spatial_rho_p0);
  }

  if (has_spatial_effects_p1 == 1) {
    REPORT(spatial_scale_p1);
    REPORT(spatial_rho_p1);
    REPORT(u_spatial_p1);
    ADREPORT(spatial_scale_p1);
    ADREPORT(spatial_rho_p1);
  }

  return nll;
}
