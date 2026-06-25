#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ============================================================================
// Kumaraswamy Distribution
// PDF: f(x) = a * b * x^(a-1) * (1 - x^a)^(b-1)
// CDF: F(x) = 1 - (1 - x^a)^b
// ============================================================================

// ----------------------------------------------------------------------------
// 1. Parametrização pela MODA (Mode)
// mu  : Moda (0 < mu < 1)
// phi : Shape 'a' (phi > 1)
//
// Conversão para 'b':
// b = (1 + (phi - 1) * mu^(-phi)) / phi
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_dkum_mode(const NumericVector x, const NumericVector mu, const NumericVector phi, const bool log_prob = false) {
  int n = std::max({(int)x.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int nx = x.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_x   = x[i % nx];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi]; // 'a'
    
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 1) { out[i] = R_NaN; continue; }
    if (cur_x <= 0 || cur_x >= 1) { out[i] = log_prob ? R_NegInf : 0.0; continue; }
    
    double cur_a = cur_phi;
    double cur_b = (1.0 + (cur_a - 1.0) * std::pow(cur_mu, -cur_a)) / cur_a;
    
    // log(f) = log(a) + log(b) + (a-1)log(x) + (b-1)log(1 - x^a)
    // log(1 - x^a) -> log1p(-x^a) para precisão
    double log_x = std::log(cur_x);
    double term_tail = std::log1p(-std::pow(cur_x, cur_a)); 
    
    double log_pdf = std::log(cur_a) + std::log(cur_b) + 
      (cur_a - 1.0) * log_x + 
      (cur_b - 1.0) * term_tail;
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pkum_mode(const NumericVector q, const NumericVector mu, const NumericVector phi, const bool lower_tail = true, const bool log_prob = false) {
  int n = std::max({(int)q.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int nq = q.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_q   = q[i % nq];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 1) { out[i] = R_NaN; continue; }
    if (cur_q <= 0) { out[i] = (log_prob ? R_NegInf : (lower_tail ? 0.0 : 1.0)); continue; }
    if (cur_q >= 1) { out[i] = (log_prob ? 0.0 : (lower_tail ? 1.0 : 0.0)); continue; }
    
    double cur_a = cur_phi;
    double cur_b = (1.0 + (cur_a - 1.0) * std::pow(cur_mu, -cur_a)) / cur_a;
    
    // F(x) = 1 - (1 - x^a)^b
    // log(1 - x^a)
    double log_base = std::log1p(-std::pow(cur_q, cur_a));
    
    // Se lower_tail = TRUE: retornar 1 - exp(b * log_base)
    // Se lower_tail = FALSE: retornar exp(b * log_base)  <-- Mais preciso calcular S(x) direto
    
    double log_S = cur_b * log_base; // log(Survival)
    
    if (lower_tail) {
      // P = 1 - S = -expm1(log_S)
      out[i] = log_prob ? std::log(-std::expm1(log_S)) : -std::expm1(log_S);
    } else {
      out[i] = log_prob ? log_S : std::exp(log_S);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qkum_mode(const NumericVector p, const NumericVector mu, const NumericVector phi, const bool lower_tail = true, const bool log_prob = false) {
  int n = std::max({(int)p.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  int np = p.size(), nmu = mu.size(), nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_p   = p[i % np];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 1) { out[i] = R_NaN; continue; }
    
    double pp = cur_p;
    if (log_prob) pp = std::exp(pp);
    if (pp < 0 || pp > 1) { out[i] = R_NaN; continue; }
    
    double cur_a = cur_phi;
    double cur_b = (1.0 + (cur_a - 1.0) * std::pow(cur_mu, -cur_a)) / cur_a;
    
    // Q(p) = (1 - (1-p)^(1/b))^(1/a)
    // Se lower_tail=FALSE: Q(p) = (1 - p^(1/b))^(1/a)
    
    double base_term; 
    if (lower_tail) {
      // (1-p)^(1/b) -> exp(1/b * log1p(-p))
      base_term = std::exp((1.0/cur_b) * std::log1p(-pp));
    } else {
      // p^(1/b) -> exp(1/b * log(p))
      base_term = std::exp((1.0/cur_b) * std::log(pp));
    }
    
    // Final: (1 - base)^(1/a) -> exp(1/a * log1p(-base))
    out[i] = std::exp((1.0/cur_a) * std::log1p(-base_term));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rkum_mode(const int n, const NumericVector mu, const NumericVector phi) {
  NumericVector u = runif(n);
  return cpp_qkum_mode(u, mu, phi, true, false);
}


// ----------------------------------------------------------------------------
// 2. Parametrização por QUANTIL (Quantile)
// mu  : Quantil de ordem tau (0 < mu < 1)
// phi : Shape 'a' (phi > 0)
// tau : Ordem do quantil (0 < tau < 1)
//
// Conversão para 'b':
// b = log(1 - tau) / log(1 - mu^a)
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector cpp_dkum_quantile(const NumericVector x, const NumericVector mu, const NumericVector phi, const NumericVector tau, const bool log_prob = false) {
  int n = std::max({(int)x.size(), (int)mu.size(), (int)phi.size(), (int)tau.size()});
  NumericVector out(n);
  int nx = x.size(), nmu = mu.size(), nphi = phi.size(), ntau = tau.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_x   = x[i % nx];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi]; // 'a'
    double cur_tau = tau[i % ntau];
    
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0 || cur_tau <= 0 || cur_tau >= 1) { out[i] = R_NaN; continue; }
    if (cur_x <= 0 || cur_x >= 1) { out[i] = log_prob ? R_NegInf : 0.0; continue; }
    
    double cur_a = cur_phi;
    // b = log(1-tau) / log(1 - mu^a)
    double cur_b = std::log1p(-cur_tau) / std::log1p(-std::pow(cur_mu, cur_a));
    
    double log_x = std::log(cur_x);
    double term_tail = std::log1p(-std::pow(cur_x, cur_a)); 
    
    double log_pdf = std::log(cur_a) + std::log(cur_b) + 
      (cur_a - 1.0) * log_x + 
      (cur_b - 1.0) * term_tail;
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pkum_quantile(const NumericVector q, const NumericVector mu, const NumericVector phi, const NumericVector tau, const bool lower_tail = true, const bool log_prob = false) {
  int n = std::max({(int)q.size(), (int)mu.size(), (int)phi.size(), (int)tau.size()});
  NumericVector out(n);
  int nq = q.size(), nmu = mu.size(), nphi = phi.size(), ntau = tau.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_q   = q[i % nq];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    double cur_tau = tau[i % ntau];
    
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0 || cur_tau <= 0 || cur_tau >= 1) { out[i] = R_NaN; continue; }
    if (cur_q <= 0) { out[i] = (log_prob ? R_NegInf : (lower_tail ? 0.0 : 1.0)); continue; }
    if (cur_q >= 1) { out[i] = (log_prob ? 0.0 : (lower_tail ? 1.0 : 0.0)); continue; }
    
    double cur_a = cur_phi;
    double cur_b = std::log1p(-cur_tau) / std::log1p(-std::pow(cur_mu, cur_a));
    
    double log_base = std::log1p(-std::pow(cur_q, cur_a));
    double log_S = cur_b * log_base;
    
    if (lower_tail) {
      out[i] = log_prob ? std::log(-std::expm1(log_S)) : -std::expm1(log_S);
    } else {
      out[i] = log_prob ? log_S : std::exp(log_S);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qkum_quantile(const NumericVector p, const NumericVector mu, const NumericVector phi, const NumericVector tau, const bool lower_tail = true, const bool log_prob = false) {
  int n = std::max({(int)p.size(), (int)mu.size(), (int)phi.size(), (int)tau.size()});
  NumericVector out(n);
  int np = p.size(), nmu = mu.size(), nphi = phi.size(), ntau = tau.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_p   = p[i % np];
    double cur_mu  = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    double cur_tau = tau[i % ntau];
    
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0 || cur_tau <= 0 || cur_tau >= 1) { out[i] = R_NaN; continue; }
    
    double pp = cur_p;
    if (log_prob) pp = std::exp(pp);
    if (pp < 0 || pp > 1) { out[i] = R_NaN; continue; }
    
    double cur_a = cur_phi;
    double cur_b = std::log1p(-cur_tau) / std::log1p(-std::pow(cur_mu, cur_a));
    
    double base_term; 
    if (lower_tail) {
      base_term = std::exp((1.0/cur_b) * std::log1p(-pp));
    } else {
      base_term = std::exp((1.0/cur_b) * std::log(pp));
    }
    out[i] = std::exp((1.0/cur_a) * std::log1p(-base_term));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rkum_quantile(const int n, const NumericVector mu, const NumericVector phi, const NumericVector tau) {
  NumericVector u = runif(n);
  return cpp_qkum_quantile(u, mu, phi, tau, true, false);
}
