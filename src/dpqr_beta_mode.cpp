#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// ============================================================================
// Beta Distribution (Mode Parameterization - Shifted Precision)
// Parameters:
//   mu  : Mode (0 < mu < 1)
//   phi : Precision (phi > 0)
//
// Relations to shapes (alpha, beta):
//   Let phi = alpha + beta - 2
//   Mode = (alpha - 1) / (alpha + beta - 2)
//   => alpha = mu * phi + 1
//   => beta  = (1 - mu) * phi + 1
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_dbeta_mode(const NumericVector x, 
                             const NumericVector mu, 
                             const NumericVector phi, 
                             const bool log_prob = false) {
  
  int n = std::max({(int)x.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int nx = x.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_x   = x[i % nx];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    // Validation: phi must be > 0 for mode to be defined in (0,1)
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) {
      out[i] = R_NaN;
      continue;
    }
    
    double shape1 = cur_mu * cur_phi + 1.0;
    double shape2 = (1.0 - cur_mu) * cur_phi + 1.0;
    
    out[i] = R::dbeta(cur_x, shape1, shape2, log_prob ? 1 : 0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pbeta_mode(const NumericVector q, 
                             const NumericVector mu, 
                             const NumericVector phi, 
                             const bool lower_tail = true, 
                             const bool log_prob = false) {
  
  int n = std::max({(int)q.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int nq = q.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_q   = q[i % nq];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) {
      out[i] = R_NaN;
      continue;
    }
    
    double shape1 = cur_mu * cur_phi + 1.0;
    double shape2 = (1.0 - cur_mu) * cur_phi + 1.0;
    
    out[i] = R::pbeta(cur_q, shape1, shape2, lower_tail ? 1 : 0, log_prob ? 1 : 0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qbeta_mode(const NumericVector p, 
                             const NumericVector mu, 
                             const NumericVector phi, 
                             const bool lower_tail = true, 
                             const bool log_prob = false) {
  
  int n = std::max({(int)p.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int np = p.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_p   = p[i % np];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) {
      out[i] = R_NaN;
      continue;
    }
    
    double shape1 = cur_mu * cur_phi + 1.0;
    double shape2 = (1.0 - cur_mu) * cur_phi + 1.0;
    
    out[i] = R::qbeta(cur_p, shape1, shape2, lower_tail ? 1 : 0, log_prob ? 1 : 0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rbeta_mode(const int n, 
                             const NumericVector mu, 
                             const NumericVector phi) {
  
  NumericVector out(n);
  int nmu = mu.size();
  int nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) {
      out[i] = R_NaN;
      continue;
    }
    
    double shape1 = cur_mu * cur_phi + 1.0;
    double shape2 = (1.0 - cur_mu) * cur_phi + 1.0;
    
    out[i] = R::rbeta(shape1, shape2);
  }
  return out;
}
