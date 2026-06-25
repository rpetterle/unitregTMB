#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ============================================================================
// Unit-Weibull Distribution (Quantile Parameterization)
// Reference: Mazucheli et al. (2020) - J. App. Stat.
//
// CDF: F(y) = tau ^ [ (log(y) / log(mu)) ^ phi ]
//
// Parameters:
//   mu  : Quantile at level tau (0 < mu < 1)
//   phi : Shape (phi > 0)
//   tau : Quantile level (0 < tau < 1)
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_duweibull(const NumericVector x, 
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
    
    // Validations
    if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0 || cur_tau <= 0 || cur_tau >= 1) { 
      out[i] = R_NaN; continue; 
    }
    if (cur_x <= 0 || cur_x >= 1) { 
      out[i] = log_prob ? R_NegInf : 0.0; continue; 
    }
    
    // Calculation using negative logs for stability
    // let ly = -log(y), lmu = -log(mu), ltau = -log(tau)
    // Ratio R = ly / lmu
    // log(f) = log(phi) + log(ltau) - log(y) + (phi-1)*log(R) - ltau * R^phi - log(lmu)
    // Simplified: log(phi) + log(ltau) - log(y) - log(lmu) + (phi-1)*(log(ly) - log(lmu)) - ltau * exp(phi * (log(ly) - log(lmu)))
    
    double nlog_y = -std::log(cur_x);     // positive
    double nlog_mu = -std::log(cur_mu);   // positive
    double nlog_tau = -std::log(cur_tau); // positive
    
    double log_ratio = std::log(nlog_y) - std::log(nlog_mu); // log(nlogy / nlogmu)
    
    double term_exp = nlog_tau * std::exp(cur_phi * log_ratio);
    
    double log_pdf = std::log(cur_phi) + std::log(nlog_tau) + nlog_y - std::log(nlog_mu) + 
      (cur_phi - 1.0) * log_ratio - term_exp;
    
    out[i] = log_prob ? log_pdf : std::exp(log_pdf);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_puweibull(const NumericVector q, 
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
    
    // F(y) = tau ^ [ (log(y)/log(mu))^phi ]
    // log(F(y)) = [ (log(y)/log(mu))^phi ] * log(tau)
    // Note: log(y)/log(mu) is positive because both are negative.
    // Let Ratio = log(y)/log(mu).
    // log(F) = - (-log(tau)) * Ratio^phi
    
    double nlog_q = -std::log(cur_q);
    double nlog_mu = -std::log(cur_mu);
    double nlog_tau = -std::log(cur_tau);
    
    double ratio = nlog_q / nlog_mu;
    // log_cdf is negative
    double log_cdf = -nlog_tau * std::pow(ratio, cur_phi);
    
    if (lower_tail) {
      out[i] = log_prob ? log_cdf : std::exp(log_cdf);
    } else {
      // S(y) = 1 - F(y) = 1 - exp(log_cdf) = -expm1(log_cdf)
      out[i] = log_prob ? std::log(-std::expm1(log_cdf)) : -std::expm1(log_cdf);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_quweibull(const NumericVector p, 
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
    // p = tau ^ [ (log(y)/log(mu))^phi ]
    // log(p) = log(tau) * [log(y)/log(mu)]^phi
    // log(p)/log(tau) = [log(y)/log(mu)]^phi
    // [log(p)/log(tau)]^(1/phi) = log(y)/log(mu)
    // log(y) = log(mu) * [log(p)/log(tau)]^(1/phi)
    // y = exp(...)
    
    double nlog_p = -std::log(pp);
    double nlog_tau = -std::log(cur_tau);
    double nlog_mu = -std::log(cur_mu);
    
    double ratio_logs = nlog_p / nlog_tau; // positive
    double term = std::pow(ratio_logs, 1.0/cur_phi);
    
    out[i] = std::exp(-nlog_mu * term);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_ruweibull(const int n, 
                            const NumericVector mu, 
                            const NumericVector phi, 
                            const NumericVector tau) {
  NumericVector u = runif(n);
  return cpp_quweibull(u, mu, phi, tau, true, false);
}