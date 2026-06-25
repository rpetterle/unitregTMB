#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ============================================================================
// Arcsecant Hyperbolic Weibull (ASHW) Distribution
// 
// Quantile Parameterization:
//   mu  : Quantile at level tau (0 < mu < 1)
//   phi : Shape parameter (phi > 0)
//   tau : Quantile level (0 < tau < 1)
//
// Auxiliary functions:
//   h(y) = arcsech(y) = log(1 + sqrt(1 - y^2)) - log(y)
//   Inverse: sech(x) = 2 / (exp(x) + exp(-x))
// ============================================================================

// Helper: arcsech(y)
// Uses difference of squares (1 - y) * (1 + y) to prevent precision loss near y = 1
inline double calc_arcsech(double y) {
  if (y <= 0.0 || y >= 1.0) return R_NaN;
  double root = std::sqrt((1.0 - y) * (1.0 + y));
  return std::log(1.0 + root) - std::log(y);
}

// Helper: sech(x)
inline double calc_sech(double x) {
  return 2.0 / (std::exp(x) + std::exp(-x));
}

// ----------------------------------------------------------------------------
// 1. PDF (Density)
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cpp_dashw(const NumericVector x, 
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
    
    double h_mu = calc_arcsech(cur_mu);
    double h_x  = calc_arcsech(cur_x);
    
    double log_h_mu = std::log(h_mu);
    double log_h_x  = std::log(h_x);
    
    // alpha = -log(tau) / (arcsech(mu)^phi)
    // Handled purely in log-space to avoid underflow/overflow
    double log_neg_log_tau = std::log(-std::log(cur_tau));
    double log_alpha = log_neg_log_tau - cur_phi * log_h_mu;
    
    // term = alpha * arcsech(x)^phi
    double log_term = log_alpha + cur_phi * log_h_x;
    double term = std::exp(log_term);
    
    // log(1 - x^2) = log(1 - x) + log(1 + x)
    double log_one_minus_x_squared = std::log(1.0 - cur_x) + std::log(1.0 + cur_x);
    
    // Final Log-PDF computation
    double log_pdf = log_alpha + std::log(cur_phi) + (cur_phi - 1.0) * log_h_x - 
                     term - std::log(cur_x) - 0.5 * log_one_minus_x_squared;
                     
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// ----------------------------------------------------------------------------
// 2. CDF (Cumulative Distribution Function)
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cpp_pashw(const NumericVector q, 
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
    
    double h_mu = calc_arcsech(cur_mu);
    double h_q  = calc_arcsech(cur_q);
    
    double log_alpha = std::log(-std::log(cur_tau)) - cur_phi * std::log(h_mu);
    double term = std::exp(log_alpha + cur_phi * std::log(h_q));
    
    // The CDF of ASHW simplifies to: F(q) = exp(-term)
    if (lower_tail) {
      out[i] = log_prob ? -term : std::exp(-term);
    } else {
      // Survival function S(q) = 1 - exp(-term) = -expm1(-term)
      out[i] = log_prob ? std::log(-std::expm1(-term)) : -std::expm1(-term);
    }
  }
  return out;
}

// ----------------------------------------------------------------------------
// 3. Quantile Function (Inverse CDF)
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cpp_qashw(const NumericVector p, 
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
    double cur_p   = p[i % np];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0 || cur_tau <= 0.0 || cur_tau >= 1.0) { 
      out[i] = R_NaN; continue; 
    }
    
    // Handle probabilities
    double pp = cur_p;
    if (log_prob) pp = std::exp(pp);
    if (!lower_tail) pp = 1.0 - pp;
    
    if (pp < 0.0 || pp > 1.0) { out[i] = R_NaN; continue; }
    if (pp == 0.0) { out[i] = 0.0; continue; }
    if (pp == 1.0) { out[i] = 1.0; continue; }
    
    // Inversion Logic:
    // p = exp(-alpha * h(y)^phi)
    // -log(p) = alpha * h(y)^phi
    // log(-log(p)) = log(alpha) + phi * log(h(y))
    // log(h(y)) = (log(-log(p)) - log(alpha)) / phi
    // Substitute log(alpha) = log(-log(tau)) - phi * log(h_mu)
    
    double log_neg_log_p   = std::log(-std::log(pp));
    double log_neg_log_tau = std::log(-std::log(cur_tau));
    double log_h_mu        = std::log(calc_arcsech(cur_mu));
    
    double log_h_y = (log_neg_log_p - log_neg_log_tau) / cur_phi + log_h_mu;
    double h_y = std::exp(log_h_y);
    
    out[i] = calc_sech(h_y);
  }
  return out;
}

// ----------------------------------------------------------------------------
// 4. Random Generation
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cpp_rashw(const int n, 
                        const NumericVector mu, 
                        const NumericVector phi, 
                        const NumericVector tau) {
  NumericVector u = runif(n);
  return cpp_qashw(u, mu, phi, tau, true, false);
}

