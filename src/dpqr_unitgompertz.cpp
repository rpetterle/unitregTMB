#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ============================================================================
// Unit-Gompertz Distribution
// Reference: Mazucheli et al. (2019)
// ============================================================================

// ----------------------------------------------------------------------------
// 1. MODE Parameterization
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_dugompertz_mode(const NumericVector x, 
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
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) { out[i] = R_NaN; continue; }
    if (cur_x <= 0.0 || cur_x >= 1.0) { out[i] = log_prob ? R_NegInf : 0.0; continue; }
    
    double alpha = std::exp(cur_phi * std::log(cur_mu)) * (cur_phi + 1.0) / cur_phi;
    
    double log_x = std::log(cur_x);
    double term = alpha * std::expm1(-cur_phi * log_x);
    double log_pdf = std::log(alpha) + std::log(cur_phi) - (cur_phi + 1.0) * log_x - term;
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pugompertz_mode(const NumericVector q, 
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
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) { out[i] = R_NaN; continue; }
    if (cur_q <= 0.0) { out[i] = (log_prob ? R_NegInf : 0.0); continue; }
    if (cur_q >= 1.0) { out[i] = (log_prob ? 0.0 : 1.0); continue; }
    
    double alpha = std::exp(cur_phi * std::log(cur_mu)) * (cur_phi + 1.0) / cur_phi;
    
    double log_cdf = -alpha * std::expm1(-cur_phi * std::log(cur_q));
    
    if (lower_tail) {
      out[i] = log_prob ? log_cdf : std::exp(log_cdf);
    } else {
      out[i] = log_prob ? std::log(-std::expm1(log_cdf)) : -std::expm1(log_cdf);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qugompertz_mode(const NumericVector p, 
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
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) { out[i] = R_NaN; continue; }
    
    double pp = p[i % np];
    if (log_prob) pp = std::exp(pp);
    if (!lower_tail) pp = 1.0 - pp;
    
    if (pp < 0.0 || pp > 1.0) { out[i] = R_NaN; continue; }
    if (pp == 0.0) { out[i] = 0.0; continue; }
    if (pp == 1.0) { out[i] = 1.0; continue; }
    
    double alpha = std::exp(cur_phi * std::log(cur_mu)) * (cur_phi + 1.0) / cur_phi;
    double term = 1.0 - (std::log(pp) / alpha);
    out[i] = std::pow(term, -1.0 / cur_phi);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rugompertz_mode(const int n, const NumericVector mu, const NumericVector phi) {
  NumericVector u = runif(n);
  return cpp_qugompertz_mode(u, mu, phi, true, false);
}

// ----------------------------------------------------------------------------
// 2. QUANTILE Parameterization
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_dugompertz_quantile(const NumericVector x, 
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
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_tau <= 0.0 || cur_tau >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    if (cur_x <= 0.0 || cur_x >= 1.0) { 
      out[i] = log_prob ? R_NegInf : 0.0; continue; 
    }
    
    double denom = std::expm1(-cur_phi * std::log(cur_mu)); 
    double alpha = -std::log(cur_tau) / denom;
    
    double log_x = std::log(cur_x);
    double term = alpha * std::expm1(-cur_phi * log_x);
    double log_pdf = std::log(alpha) + std::log(cur_phi) - (cur_phi + 1.0) * log_x - term;
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pugompertz_quantile(const NumericVector q, 
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
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_tau <= 0.0 || cur_tau >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    if (cur_q <= 0.0) { out[i] = (log_prob ? R_NegInf : 0.0); continue; }
    if (cur_q >= 1.0) { out[i] = (log_prob ? 0.0 : 1.0); continue; }
    
    double denom = std::expm1(-cur_phi * std::log(cur_mu)); 
    double alpha = -std::log(cur_tau) / denom;
    
    double log_cdf = -alpha * std::expm1(-cur_phi * std::log(cur_q));
    
    if (lower_tail) {
      out[i] = log_prob ? log_cdf : std::exp(log_cdf);
    } else {
      out[i] = log_prob ? std::log(-std::expm1(log_cdf)) : -std::expm1(log_cdf);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qugompertz_quantile(const NumericVector p, 
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
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_tau <= 0.0 || cur_tau >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    
    double pp = p[i % np];
    if (log_prob) pp = std::exp(pp);
    if (!lower_tail) pp = 1.0 - pp;
    
    if (pp < 0.0 || pp > 1.0) { out[i] = R_NaN; continue; }
    if (pp == 0.0) { out[i] = 0.0; continue; }
    if (pp == 1.0) { out[i] = 1.0; continue; }
    
    double log_tau = std::log(cur_tau);
    double denom = std::expm1(-cur_phi * std::log(cur_mu)); 
    
    double term = 1.0 + (std::log(pp) / log_tau) * denom;
    out[i] = std::pow(term, -1.0 / cur_phi);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rugompertz_quantile(const int n, 
                                      const NumericVector mu, 
                                      const NumericVector phi, 
                                      const NumericVector tau) {
  NumericVector u = runif(n);
  return cpp_qugompertz_quantile(u, mu, phi, tau, true, false);
}
