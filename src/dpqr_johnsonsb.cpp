#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ============================================================================
// Johnson SB Distribution (Quantile Parameterization)
// Transformation: Z = gamma + delta * log(Y / (1-Y)) ~ N(0,1)
//
// Reparameterized in terms of Quantile (mu) at level (tau):
//   Z = z_tau + phi * (logit(y) - logit(mu))
// where z_tau = qnorm(tau)
//
// Parameters:
//   mu  : Quantile parameter (0 < mu < 1)
//   phi : Shape parameter (phi > 0)
//   tau : Quantile level (0 < tau < 1)
// ============================================================================

// Helper: Logit function
inline double logit(double p) {
  return std::log(p / (1.0 - p));
}

// Helper: Inverse Logit (Expit)
inline double expit(double z) {
  return 1.0 / (1.0 + std::exp(-z));
}

// [[Rcpp::export]]
NumericVector cpp_djohnsonsb(const NumericVector x, 
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
    
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0 || cur_tau <= 0 || cur_tau >= 1) { 
      out[i] = R_NaN; continue; 
    }
    if (cur_x <= 0 || cur_x >= 1) { 
      out[i] = log_prob ? R_NegInf : 0.0; continue; 
    }
    
    // Z-score calculation
    double z_tau = R::qnorm(cur_tau, 0.0, 1.0, 1, 0);
    double logit_x = logit(cur_x);
    double logit_mu = logit(cur_mu);
    
    // Linear predictor in normal scale
    double z = z_tau + cur_phi * (logit_x - logit_mu);
    
    // Density of standard normal at z
    double dnorm_z = R::dnorm(z, 0.0, 1.0, 1); // log scale
    
    // Jacobian of transformation: phi / (x * (1-x))
    // log(Jacobian) = log(phi) - log(x) - log(1-x)
    double log_jac = std::log(cur_phi) - std::log(cur_x) - std::log(1.0 - cur_x);
    
    double log_pdf = dnorm_z + log_jac;
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pjohnsonsb(const NumericVector q, 
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
    
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0 || cur_tau <= 0 || cur_tau >= 1) { 
      out[i] = R_NaN; continue; 
    }
    if (cur_q <= 0) { out[i] = (log_prob ? R_NegInf : 0.0); continue; }
    if (cur_q >= 1) { out[i] = (log_prob ? 0.0 : 1.0); continue; }
    
    double z_tau = R::qnorm(cur_tau, 0.0, 1.0, 1, 0);
    double logit_q = logit(cur_q);
    double logit_mu = logit(cur_mu);
    
    double z = z_tau + cur_phi * (logit_q - logit_mu);
    
    out[i] = R::pnorm(z, 0.0, 1.0, lower_tail ? 1 : 0, log_prob ? 1 : 0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qjohnsonsb(const NumericVector p, 
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
    
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0 || cur_tau <= 0 || cur_tau >= 1) { 
      out[i] = R_NaN; continue; 
    }
    
    // Handle probabilities
    double pp = cur_p;
    if (log_prob) pp = std::exp(pp);
    if (!lower_tail) pp = 1.0 - pp;
    
    if (pp < 0 || pp > 1) { out[i] = R_NaN; continue; }
    if (pp == 0) { out[i] = 0.0; continue; }
    if (pp == 1) { out[i] = 1.0; continue; }
    
    // Inversion logic
    // p = Phi(z_tau + phi * (logit(y) - logit(mu)))
    // z_p = z_tau + phi * (logit(y) - logit(mu))
    // (z_p - z_tau) / phi = logit(y) - logit(mu)
    // logit(y) = logit(mu) + (z_p - z_tau) / phi
    // y = expit(...)
    
    double z_tau = R::qnorm(cur_tau, 0.0, 1.0, 1, 0);
    double z_p   = R::qnorm(pp, 0.0, 1.0, 1, 0);
    double logit_mu = logit(cur_mu);
    
    double target_logit = logit_mu + (z_p - z_tau) / cur_phi;
    
    out[i] = expit(target_logit);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rjohnsonsb(const int n, 
                             const NumericVector mu, 
                             const NumericVector phi, 
                             const NumericVector tau) {
  NumericVector u = runif(n);
  return cpp_qjohnsonsb(u, mu, phi, tau, true, false);
}

