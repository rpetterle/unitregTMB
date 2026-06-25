#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ============================================================================
// Unit-Birnbaum-Saunders (UBS) Distribution
// Reference: Mazucheli et al. (2020)
//
// Transformation: Y = exp(-T), T ~ BS(alpha, beta)
//
// Parameters:
//   mu  : Quantile at level tau (0 < mu < 1)
//   phi : Shape (phi > 0)
//   tau : Quantile level (0 < tau < 1)
// ============================================================================

// Helper: t(y) = sqrt(-1 / log(y))
inline double calc_t(double y) {
  if (y <= 0 || y >= 1) return R_NaN;
  return std::sqrt(-1.0 / std::log(y));
}

// Helper: xi(y) = t(y) - 1/t(y)
inline double calc_xi(double t_val) {
  return t_val - 1.0 / t_val;
}

// [[Rcpp::export]]
NumericVector cpp_dubs(const NumericVector x, 
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
    
    // 1. Calculate auxiliary terms
    double t_x = calc_t(cur_x);
    double t_mu = calc_t(cur_mu);
    
    double xi_x = calc_xi(t_x);
    double xi_mu = calc_xi(t_mu);
    
    double z_tau = R::qnorm(cur_tau, 0.0, 1.0, 1, 0); // lower_tail=1, log=0
    
    // 2. Z-score
    double z = (xi_x - xi_mu) / cur_phi + z_tau;
    
    // 3. Jacobian components
    // J = (t(x) + 1/t(x)) / (2 * x * (-log(x)) * phi)
    double term_num = t_x + 1.0 / t_x;
    double term_den = 2.0 * cur_x * (-std::log(cur_x)) * cur_phi;
    
    // 4. Log-PDF
    // log(f) = dnorm(z, log=T) + log(J)
    double log_dnorm = R::dnorm(z, 0.0, 1.0, 1);
    double log_jac = std::log(term_num) - std::log(term_den);
    
    double log_pdf = log_dnorm + log_jac;
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pubs(const NumericVector q, 
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
    
    // F(y) = Phi(z)
    double t_q = calc_t(cur_q);
    double t_mu = calc_t(cur_mu);
    
    double xi_q = calc_xi(t_q);
    double xi_mu = calc_xi(t_mu);
    double z_tau = R::qnorm(cur_tau, 0.0, 1.0, 1, 0);
    
    double z = (xi_q - xi_mu) / cur_phi + z_tau;
    
    out[i] = R::pnorm(z, 0.0, 1.0, lower_tail ? 1 : 0, log_prob ? 1 : 0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qubs(const NumericVector p, 
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
    
    double pp = cur_p;
    if (log_prob) pp = std::exp(pp);
    if (!lower_tail) pp = 1.0 - pp;
    
    if (pp < 0 || pp > 1) { out[i] = R_NaN; continue; }
    if (pp == 0) { out[i] = 0.0; continue; }
    if (pp == 1) { out[i] = 1.0; continue; }
    
    // Inversion:
    // z_p = Phi^-1(p)
    // z_p = (xi(y) - xi(mu))/phi + z_tau
    // xi(y) = phi * (z_p - z_tau) + xi(mu)
    // Let K = phi * (z_p - z_tau) + xi(mu)
    // xi(y) = t(y) - 1/t(y) = K
    // t^2 - K*t - 1 = 0
    // t = (K + sqrt(K^2 + 4)) / 2  (Only positive root valid)
    // log(y) = -1 / t^2
    // y = exp(-1 / t^2)
    
    double z_p = R::qnorm(pp, 0.0, 1.0, 1, 0);
    double z_tau = R::qnorm(cur_tau, 0.0, 1.0, 1, 0);
    
    double t_mu = calc_t(cur_mu);
    double xi_mu = calc_xi(t_mu);
    
    double K = cur_phi * (z_p - z_tau) + xi_mu;
    
    double t_y = (K + std::sqrt(K*K + 4.0)) / 2.0;
    
    out[i] = std::exp(-1.0 / (t_y * t_y));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rubs(const int n, 
                       const NumericVector mu, 
                       const NumericVector phi, 
                       const NumericVector tau) {
  NumericVector u = runif(n);
  return cpp_qubs(u, mu, phi, tau, true, false);
}