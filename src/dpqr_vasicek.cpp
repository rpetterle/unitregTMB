#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ============================================================================
// Vasicek Distribution
// Reference: Mazucheli et al. (2022) Parameterizations
// 
// Important Note: In this package architecture, the parameter 'phi' 
// represents the standard shape parameter 'theta' bounded in (0, 1). 
// The logit link transformation is handled upstream by the TMB model environment.
// ============================================================================

// ----------------------------------------------------------------------------
// 1. MEAN Parameterization
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_dvasicek_mean(const NumericVector x, 
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
    
    // Bounds check: phi is strictly bounded in (0, 1)
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_phi >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    if (cur_x <= 0.0 || cur_x >= 1.0) { 
      out[i] = log_prob ? R_NegInf : 0.0; continue; 
    }
    
    double z_y  = R::qnorm(cur_x, 0.0, 1.0, 1, 0);
    double z_mu = R::qnorm(cur_mu, 0.0, 1.0, 1, 0);
    
    double arg = (std::sqrt(1.0 - cur_phi) * z_y - z_mu) / std::sqrt(cur_phi);
    double log_pdf = 0.5 * (std::log(1.0 - cur_phi) - std::log(cur_phi) + z_y * z_y - arg * arg);
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pvasicek_mean(const NumericVector q, 
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
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_phi >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    if (cur_q <= 0.0) { out[i] = (log_prob ? R_NegInf : 0.0); continue; }
    if (cur_q >= 1.0) { out[i] = (log_prob ? 0.0 : 1.0); continue; }
    
    double z_q  = R::qnorm(cur_q, 0.0, 1.0, 1, 0);
    double z_mu = R::qnorm(cur_mu, 0.0, 1.0, 1, 0);
    
    double arg = (std::sqrt(1.0 - cur_phi) * z_q - z_mu) / std::sqrt(cur_phi);
    out[i] = R::pnorm(arg, 0.0, 1.0, lower_tail ? 1 : 0, log_prob ? 1 : 0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qvasicek_mean(const NumericVector p, 
                                const NumericVector mu, 
                                const NumericVector phi, 
                                const bool lower_tail = true, 
                                const bool log_prob = false) {
  int n = std::max({(int)p.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int np = p.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_phi >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    
    double cur_p = p[i % np];
    double prob = log_prob ? std::exp(cur_p) : cur_p;
    if (!lower_tail) prob = 1.0 - prob;
    
    if (prob < 0.0 || prob > 1.0) { out[i] = R_NaN; continue; }
    if (prob == 0.0) { out[i] = 0.0; continue; }
    if (prob == 1.0) { out[i] = 1.0; continue; }
    
    double z_p  = R::qnorm(cur_p, 0.0, 1.0, lower_tail ? 1 : 0, log_prob ? 1 : 0);
    double z_mu = R::qnorm(cur_mu, 0.0, 1.0, 1, 0);
    
    double arg = (z_mu + std::sqrt(cur_phi) * z_p) / std::sqrt(1.0 - cur_phi);
    out[i] = R::pnorm(arg, 0.0, 1.0, 1, 0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rvasicek_mean(const int n, 
                                const NumericVector mu, 
                                const NumericVector phi) {
  NumericVector u = runif(n);
  return cpp_qvasicek_mean(u, mu, phi, true, false);
}


// ----------------------------------------------------------------------------
// 2. QUANTILE Parameterization
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_dvasicek_quantile(const NumericVector x, 
                                    const NumericVector mu, 
                                    const NumericVector phi, 
                                    const NumericVector tau, 
                                    const bool log_prob = false) {
  int n = std::max({(int)x.size(), (int)mu.size(), (int)phi.size(), (int)tau.size()});
  NumericVector out(n);
  int nx = x.size(), nmu = mu.size(), nphi = phi.size(), ntau = tau.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_x   = x[i % nx];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    double cur_tau = tau[i % ntau];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_phi >= 1.0 || cur_tau <= 0.0 || cur_tau >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    if (cur_x <= 0.0 || cur_x >= 1.0) { 
      out[i] = log_prob ? R_NegInf : 0.0; continue; 
    }
    
    double z_y   = R::qnorm(cur_x, 0.0, 1.0, 1, 0);
    double z_mu  = R::qnorm(cur_mu, 0.0, 1.0, 1, 0);
    double z_tau = R::qnorm(cur_tau, 0.0, 1.0, 1, 0);
    
    double num = std::sqrt(1.0 - cur_phi) * (z_y - z_mu) + std::sqrt(cur_phi) * z_tau;
    double arg = num / std::sqrt(cur_phi);
    
    double log_pdf = 0.5 * (std::log(1.0 - cur_phi) - std::log(cur_phi) + z_y * z_y - arg * arg);
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pvasicek_quantile(const NumericVector q, 
                                    const NumericVector mu, 
                                    const NumericVector phi, 
                                    const NumericVector tau, 
                                    const bool lower_tail = true, 
                                    const bool log_prob = false) {
  int n = std::max({(int)q.size(), (int)mu.size(), (int)phi.size(), (int)tau.size()});
  NumericVector out(n);
  int nq = q.size(), nmu = mu.size(), nphi = phi.size(), ntau = tau.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_q   = q[i % nq];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    double cur_tau = tau[i % ntau];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_phi >= 1.0 || cur_tau <= 0.0 || cur_tau >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    if (cur_q <= 0.0) { out[i] = (log_prob ? R_NegInf : 0.0); continue; }
    if (cur_q >= 1.0) { out[i] = (log_prob ? 0.0 : 1.0); continue; }
    
    double z_q   = R::qnorm(cur_q, 0.0, 1.0, 1, 0);
    double z_mu  = R::qnorm(cur_mu, 0.0, 1.0, 1, 0);
    double z_tau = R::qnorm(cur_tau, 0.0, 1.0, 1, 0);
    
    double arg = std::sqrt((1.0 - cur_phi) / cur_phi) * (z_q - z_mu) + z_tau;
    out[i] = R::pnorm(arg, 0.0, 1.0, lower_tail ? 1 : 0, log_prob ? 1 : 0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qvasicek_quantile(const NumericVector p, 
                                    const NumericVector mu, 
                                    const NumericVector phi, 
                                    const NumericVector tau, 
                                    const bool lower_tail = true, 
                                    const bool log_prob = false) {
  int n = std::max({(int)p.size(), (int)mu.size(), (int)phi.size(), (int)tau.size()});
  NumericVector out(n);
  int np = p.size(), nmu = mu.size(), nphi = phi.size(), ntau = tau.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    double cur_tau = tau[i % ntau];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_phi >= 1.0 || cur_tau <= 0.0 || cur_tau >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    
    double cur_p = p[i % np];
    double prob = log_prob ? std::exp(cur_p) : cur_p;
    if (!lower_tail) prob = 1.0 - prob;
    
    if (prob < 0.0 || prob > 1.0) { out[i] = R_NaN; continue; }
    if (prob == 0.0) { out[i] = 0.0; continue; }
    if (prob == 1.0) { out[i] = 1.0; continue; }
    
    double z_p   = R::qnorm(cur_p, 0.0, 1.0, lower_tail ? 1 : 0, log_prob ? 1 : 0);
    double z_mu  = R::qnorm(cur_mu, 0.0, 1.0, 1, 0);
    double z_tau = R::qnorm(cur_tau, 0.0, 1.0, 1, 0);
    
    double arg = z_mu + std::sqrt(cur_phi / (1.0 - cur_phi)) * (z_p - z_tau);
    out[i] = R::pnorm(arg, 0.0, 1.0, 1, 0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rvasicek_quantile(const int n, 
                                    const NumericVector mu, 
                                    const NumericVector phi, 
                                    const NumericVector tau) {
  NumericVector u = runif(n);
  return cpp_qvasicek_quantile(u, mu, phi, tau, true, false);
}

