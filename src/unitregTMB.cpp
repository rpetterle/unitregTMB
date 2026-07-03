// ----------------------------------------------------------------------------
// unitregTMB.cpp 
// ----------------------------------------------------------------------------

#define EIGEN_DONT_PARALLELIZE 
#define TMB_LIB_INIT tmb_custom_init

#include <TMB.hpp>
#include <cmath> 

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "Distributions.hpp" 


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
    case logit_link:   ans = invlogit(eta); break;
    case probit_link:  ans = pnorm(eta); break;
    case cloglog_link: ans = Type(1.0) - exp(-exp(eta)); break;
    case cauchit_link: ans = Type(0.5) + atan(eta) / Type(M_PI); break;
    default: error("Link not implemented or invalid for mu!");
  }
  return ans;
}

template<class Type>
Type get_log_density(int family, Type y, Type mu, Type phi, Type tau) {
  switch (family) {
    // Mean (4 arguments) 
    case fam_beta_mean:          return dbeta_mean(y, mu, phi, 1);
    case fam_simplex_mean:       return dsimplex(y, mu, phi, 1); 
    case fam_vasicek_mean:       return dvasicek_mean(y, mu, phi, 1);
    case fam_ugamma_mean:        return dunitgamma_mean(y, mu, phi, 1);
    case fam_bessel_mean:        return dbessel_mean(y, mu, phi, 1);
    
    // Mode (4 arguments) 
    case fam_beta_mode:          return dbeta_mode(y, mu, phi, 1);
    case fam_kum_mode:           return dkum_mode(y, mu, phi, 1);
    case fam_ugamma_mode:        return dugamma_mode(y, mu, phi, 1);
    case fam_ugompertz_mode:     return dunitgompertz_mode(y, mu, phi, 1);
    
    // Quantile (5 arguments: include tau) 
    case fam_kum_quantile:       return dkum_quantile(y, mu, phi, tau, 1);
    case fam_vasicek_quantile:   return dvasicek_quantile(y, mu, phi, tau, 1);
    case fam_uweibull_quantile:  return dunitweibull(y, mu, phi, tau, 1);
    case fam_ugompertz_quantile: return dunitgompertz_quantile(y, mu, phi, tau, 1);
    case fam_johnsonsb_quantile: return djohnsonsb(y, mu, phi, tau, 1);
    case fam_ashw_quantile:      return dashw(y, mu, phi, tau, 1);
    case fam_ubs_quantile:       return dubs(y, mu, phi, tau, 1);
    
    default: error("Unsupported family code."); return Type(0.0);
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


template<class Type>
Type objective_function<Type>::operator()() {
  
  DATA_VECTOR(Y);
  DATA_IVECTOR(y_class);   
  DATA_MATRIX(X_mu);                        
  DATA_MATRIX(X_phi);                       
  DATA_MATRIX(X_p0);                        
  DATA_MATRIX(X_p1);                        
  DATA_SPARSE_MATRIX(Z_mu);
  DATA_VECTOR(weights);                     
  DATA_VECTOR(offset_mu);                 
  DATA_INTEGER(family);                     
  DATA_INTEGER(link);                       
  DATA_INTEGER(has_p0_inflation);           
  DATA_INTEGER(has_p1_inflation);           
  DATA_INTEGER(has_random_effects_mu);      
  DATA_SCALAR(tau);                         

  DATA_IVECTOR(re_term_n_components);       
  DATA_IVECTOR(re_term_chol_starts);        
  DATA_IVECTOR(u_mu_term_starts);           
  DATA_IVECTOR(n_re_levels_list_data);      
  
  PARAMETER_VECTOR(beta_mu);                
  PARAMETER_VECTOR(beta_phi);               
  PARAMETER_VECTOR(beta_p0);                
  PARAMETER_VECTOR(beta_p1);                
  PARAMETER_VECTOR(u_mu);                   
  PARAMETER_VECTOR(log_chol_re_mu_combined);
  
  Type nll = Type(0.0);
  int n_obs = Y.size();
  
  vector<Type> conditional_log_lik_unweighted(n_obs);
  vector<Type> conditional_log_lik_weighted(n_obs);
  conditional_log_lik_unweighted.setZero();
  conditional_log_lik_weighted.setZero();
  
  vector<Type> eta_mu = X_mu * beta_mu + offset_mu;
  if (has_random_effects_mu == 1) {
    eta_mu += Z_mu * u_mu;
  }
  
  vector<Type> eta_phi = X_phi * beta_phi;
  
  vector<Type> eta_p0(n_obs);
  eta_p0.setZero();
  if (has_p0_inflation == 1) {
    eta_p0 = X_p0 * beta_p0;
  } 
  
  vector<Type> eta_p1(n_obs);
  eta_p1.setZero();
  if (has_p1_inflation == 1) {
    eta_p1 = X_p1 * beta_p1;
  } 

  int total_n_comp = 0;
  int total_n_cor = 0;
  if (has_random_effects_mu == 1) {
    for (int t = 0; t < re_term_n_components.size(); ++t) {
        int nc = re_term_n_components(t);
        total_n_comp += nc;
        total_n_cor += nc * (nc - 1) / 2;
    }
  }
  
  vector<Type> re_std_dev(total_n_comp);
  vector<Type> re_cor(total_n_cor > 0 ? total_n_cor : 1); 
  int idx_std = 0;
  int idx_cor = 0;

  if (has_random_effects_mu == 1) {
    using namespace density;
    int n_re_terms = re_term_n_components.size();
    
    for (int t_idx = 0; t_idx < n_re_terms; ++t_idx) {
      int n_comp = re_term_n_components(t_idx);
      int chol_start = re_term_chol_starts(t_idx);
      int u_start = u_mu_term_starts(t_idx);
      int n_levels = n_re_levels_list_data(t_idx);
      
      matrix<Type> L(n_comp, n_comp);
      L.setZero();
      int k_param = 0;
      for(int j=0; j<n_comp; ++j){
        for(int i=j; i<n_comp; ++i){
          if(i==j) L(i,j) = exp(log_chol_re_mu_combined(chol_start + k_param));
          else     L(i,j) = log_chol_re_mu_combined(chol_start + k_param);
          k_param++;
        }
      }
      
      matrix<Type> Sigma = L * L.transpose();
      for(int j = 0; j < n_comp; ++j) {
        re_std_dev(idx_std++) = sqrt(Sigma(j,j));
        for(int i = j + 1; i < n_comp; ++i) {
           re_cor(idx_cor++) = Sigma(i,j) / sqrt(Sigma(i,i) * Sigma(j,j));
        }
      }
      
      MVNORM_t<Type> neg_log_dmvnorm(Sigma);
      
      for(int k = 0; k < n_levels; ++k) {
        vector<Type> u_block = u_mu.segment(u_start + k * n_comp, n_comp);
        nll += neg_log_dmvnorm(u_block); 
      }
    }
  }
  
  for(int i = 0; i < n_obs; ++i) {
      
      Type mu_i  = inverse_linkfun(eta_mu(i), link); 
      
      Type phi_i = inverse_phi_linkfun(eta_phi(i), family);

      Type log_p0 = Type(0.0);
      Type log_p1 = Type(0.0);
      Type log_pc = Type(0.0);
      Type li = Type(0.0); 
      
      if (has_p0_inflation == 1 && has_p1_inflation == 1) {
          Type log_boundary_sum = logspace_add(eta_p0(i), eta_p1(i));
          Type log_denom        = logspace_add(Type(0.0), log_boundary_sum);
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
  
  REPORT(conditional_log_lik_unweighted);
  REPORT(conditional_log_lik_weighted);
  
  if (has_random_effects_mu == 1) {
    REPORT(re_std_dev);        
    if (total_n_cor > 0) {
      REPORT(re_cor);          
    }
    
    ADREPORT(log_chol_re_mu_combined);
    ADREPORT(re_std_dev);        
    if (total_n_cor > 0) {
      ADREPORT(re_cor);          
    }
  }
  
  return nll;
}
