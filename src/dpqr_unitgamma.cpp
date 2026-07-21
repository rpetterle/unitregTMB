#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ============================================================================
// Unit-Gamma Distribution
// Relation: Y = exp(-X), where X ~ Gamma(shape=phi, rate=beta)
// ============================================================================

// ----------------------------------------------------------------------------
// 1. MEAN Parameterization
// mu  : Mean (0 < mu < 1)
// phi : Shape (phi > 0)
//
// Conversion to Rate (beta):
// beta = 1 / (mu^(-1/phi) - 1)
// Using expm1 for numerical stability when mu approaches 1.
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_dugamma_mean(const NumericVector x, const NumericVector mu, const NumericVector phi, const bool log_prob = false) {
  int n = std::max({(int)x.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int nx = x.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_x   = x[i % nx];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) { out[i] = R_NaN; continue; }
    if (cur_x <= 0.0 || cur_x >= 1.0) { out[i] = log_prob ? R_NegInf : 0.0; continue; }
    
    double a = -std::log(cur_mu) / cur_phi;
    double beta = 1.0 / std::expm1(a);
    
    double log_y = std::log(cur_x);
    double log_neg_log_y = std::log(-log_y);
    
    // Density: f(y) = (beta^phi / Gamma(phi)) * y^(beta-1) * (-log(y))^(phi-1)
    double log_pdf = cur_phi * std::log(beta) - R::lgammafn(cur_phi) + 
                     (beta - 1.0) * log_y + (cur_phi - 1.0) * log_neg_log_y;
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pugamma_mean(const NumericVector q, const NumericVector mu, const NumericVector phi, const bool lower_tail = true, const bool log_prob = false) {
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
    
    double a = -std::log(cur_mu) / cur_phi;
    double beta = 1.0 / std::expm1(a);
    
    out[i] = R::pgamma(-std::log(cur_q), cur_phi, 1.0/beta, !lower_tail, log_prob);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qugamma_mean(const NumericVector p, const NumericVector mu, const NumericVector phi, const bool lower_tail = true, const bool log_prob = false) {
  int n = std::max({(int)p.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int np = p.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    double cur_p   = p[i % np];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) { out[i] = R_NaN; continue; }
    
    double pp = cur_p;
    if (log_prob) pp = std::exp(pp);
    if (pp < 0.0 || pp > 1.0) { out[i] = R_NaN; continue; }
    if (pp == 0.0) { out[i] = lower_tail ? 0.0 : 1.0; continue; }
    if (pp == 1.0) { out[i] = lower_tail ? 1.0 : 0.0; continue; }
    
    double a = -std::log(cur_mu) / cur_phi;
    double beta = 1.0 / std::expm1(a);
    
    double x_quant = R::qgamma(pp, cur_phi, 1.0/beta, !lower_tail, false);
    out[i] = std::exp(-x_quant);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rugamma_mean(const int n, const NumericVector mu, const NumericVector phi) {
  NumericVector out(n);
  int nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) { out[i] = R_NaN; continue; }
    
    double a = -std::log(cur_mu) / cur_phi;
    double beta = 1.0 / std::expm1(a);
    
    out[i] = std::exp(-R::rgamma(cur_phi, 1.0/beta));
  }
  return out;
}


// ----------------------------------------------------------------------------
// 2. MODE Parameterization
// mu  : Mode (0 < mu < 1)
// phi : Shape (phi > 0)
//
// Conversion to Rate (beta):
// beta = (1 + log(mu) - phi) / log(mu)
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_dugamma_mode(const NumericVector x, const NumericVector mu, const NumericVector phi, const bool log_prob = false) {
  int n = std::max({(int)x.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int nx = x.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_x   = x[i % nx];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) { out[i] = R_NaN; continue; }
    if (cur_x <= 0.0 || cur_x >= 1.0) { out[i] = log_prob ? R_NegInf : 0.0; continue; }
    
    double log_mu = std::log(cur_mu);
    double beta = (1.0 + log_mu - cur_phi) / log_mu;
    
    if (beta <= 0.0) { out[i] = R_NaN; continue; }
    
    double log_y = std::log(cur_x);
    double log_neg_log_y = std::log(-log_y);
    
    double log_pdf = cur_phi * std::log(beta) - R::lgammafn(cur_phi) + 
                     (beta - 1.0) * log_y + (cur_phi - 1.0) * log_neg_log_y;
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pugamma_mode(const NumericVector q, const NumericVector mu, const NumericVector phi, const bool lower_tail = true, const bool log_prob = false) {
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
    
    double log_mu = std::log(cur_mu);
    double beta = (1.0 + log_mu - cur_phi) / log_mu;
    if (beta <= 0.0) { out[i] = R_NaN; continue; }
    
    out[i] = R::pgamma(-std::log(cur_q), cur_phi, 1.0/beta, !lower_tail, log_prob);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qugamma_mode(const NumericVector p, const NumericVector mu, const NumericVector phi, const bool lower_tail = true, const bool log_prob = false) {
  int n = std::max({(int)p.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int np = p.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    double cur_p   = p[i % np];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) { out[i] = R_NaN; continue; }
    
    double pp = cur_p;
    if (log_prob) pp = std::exp(pp);
    if (pp < 0.0 || pp > 1.0) { out[i] = R_NaN; continue; }
    if (pp == 0.0) { out[i] = lower_tail ? 0.0 : 1.0; continue; }
    if (pp == 1.0) { out[i] = lower_tail ? 1.0 : 0.0; continue; }
    
    double log_mu = std::log(cur_mu);
    double beta = (1.0 + log_mu - cur_phi) / log_mu;
    if (beta <= 0.0) { out[i] = R_NaN; continue; }
    
    double x_quant = R::qgamma(pp, cur_phi, 1.0/beta, !lower_tail, false);
    out[i] = std::exp(-x_quant);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rugamma_mode(const int n, const NumericVector mu, const NumericVector phi) {
  NumericVector out(n);
  int nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) { out[i] = R_NaN; continue; }
    
    double log_mu = std::log(cur_mu);
    double beta = (1.0 + log_mu - cur_phi) / log_mu;
    if (beta <= 0.0) { out[i] = R_NaN; continue; }
    
    out[i] = std::exp(-R::rgamma(cur_phi, 1.0/beta));
  }
  return out;
}
