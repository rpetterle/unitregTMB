#ifndef UNITREG_DISTRIBUTIONS_HPP
#define UNITREG_DISTRIBUTIONS_HPP

// Helper functions -----------------------------------------------------------
template<class Type>
Type arcsech(Type y) {
  Type one = Type(1.0);
  Type root = sqrt((one - y) * (one + y));
  return log(one + root) - log(y);
}

// Parametrization: mean ------------------------------------------------------ 

// 0. Beta (mean)
template<class Type>
Type dbeta_mean(Type y, Type mu, Type phi, int give_log) {
  Type shape1 = mu * phi;
  Type shape2 = (Type(1.0) - mu) * phi;
  return dbeta(y, shape1, shape2, give_log);
}

// 1. Simplex (mean)
template<class Type>
Type dsimplex(Type y, Type mu, Type phi, int give_log) {
  Type one = Type(1.0);
  Type half = Type(0.5);
  Type one_minus_y  = one - y;
  Type one_minus_mu = one - mu;
  Type diff = y - mu;
  Type d = (diff * diff) / (y * one_minus_y * mu * mu * one_minus_mu * one_minus_mu);
  Type log_res = -half * (log(Type(2.0) * Type(M_PI)) + log(phi) + 
                Type(3.0) * (log(y) + log(one_minus_y))) - d / (Type(2.0) * phi);
  
  return give_log ? log_res : exp(log_res);
}

// 2. Vasicek (mean)
template<class Type>
Type dvasicek_mean(Type y, Type mu, Type phi, int give_log) {
  Type one  = Type(1.0);
  Type half = Type(0.5);
  Type one_minus_phi = one - phi;
  
  Type z_y  = qnorm(y);
  Type z_mu = qnorm(mu);
  
  Type arg = (sqrt(one_minus_phi) * z_y - z_mu) / sqrt(phi);
  Type log_pdf = half * (log(one_minus_phi) - log(phi) + z_y * z_y - arg * arg);
  
  return give_log ? log_pdf : exp(log_pdf);
}

// 3. Unit-gamma (mean)
template<class Type>
Type dunitgamma_mean(Type y, Type mu, Type phi, int give_log) {
  Type one  = Type(1.0);
  Type zero = Type(0.0);
  
  Type log_y = log(y);
  Type log_mu = log(mu);
  Type nlog_y = -log_y;
  Type a = -log_mu / phi;
  Type log_denom = logspace_sub(a, zero);
  Type log_beta = -log_denom;
  Type beta = exp(log_beta);
  Type log_pdf = phi * log_beta - lgamma(phi) + (beta - one) * log_y +
                 (phi - one) * log(nlog_y);
  
  return give_log ? log_pdf : exp(log_pdf);
}

// 4. Bessel (mean)
template<class Type>
Type dbessel_mean(Type y, Type mu, Type phi, int give_log) {
  Type one  = Type(1.0);
  Type half = Type(0.5);
  Type eps  = Type(1e-12);
  Type y_safe = (y < eps) ? eps : ((y > one - eps) ? (one - eps) : y);
  Type log_y = log(y_safe);
  Type log_one_minus_y = log(one - y_safe);
  Type den_zeta = y_safe * (one - y_safe);
  Type diff = y_safe - mu;
  Type zeta = sqrt(one + (diff * diff) / den_zeta);
  Type arg = phi * zeta;
  Type log_k1 = log(besselK(arg, one));
  
  Type log_num = log(mu) + log(one - mu) + log(phi) + phi;
  Type log_den = log(Type(M_PI)) + Type(1.5) * (log_y + log_one_minus_y) + log(zeta);
  Type log_res = log_num + log_k1 - log_den;
 
  return give_log ? log_res : exp(log_res);
}

// 5. Beta (mode)
template<class Type>
Type dbeta_mode(Type y, Type mu, Type phi, int give_log) {
  Type shape1 = mu * phi + Type(1.0);
  Type shape2 = (Type(1.0) - mu) * phi + Type(1.0);
  return dbeta(y, shape1, shape2, give_log);
}

// 6. Kumaraswamy (mode)
template<class Type>
Type dkum_mode(Type x, Type mu, Type phi, int give_log) {
  Type one = Type(1.0);
  Type log_x  = log(x);
  Type log_mu = log(mu);
  Type mu_inv_phi = exp(-phi * log_mu);
  Type b = (one + (phi - one) * mu_inv_phi) / phi;
  Type log_x_phi = phi * log_x;
  Type log_one_minus_xphi = logspace_sub(Type(0.0), log_x_phi);
  Type log_pdf = log(phi) + log(b) + (phi - one) * log_x +
    (b - one) * log_one_minus_xphi;
  return give_log ? log_pdf : exp(log_pdf);
}

// 7. Unit-gamma (mode)
template<class Type>
Type dugamma_mode(Type y, Type mu, Type phi, int give_log) {
  Type one = Type(1.0);
  Type log_y = log(y);
  Type nlog_y = -log_y;
  Type nlog_mu = -log(mu);
  Type beta = one + (phi - one) / nlog_mu;
  
  Type log_pdf = phi * log(beta) - lgamma(phi) + (beta - one) * log_y +
                 (phi - one) * log(nlog_y);
  return give_log ? log_pdf : exp(log_pdf);
}

// 8. Unit-Gompertz (mode)
template<class Type>
Type dunitgompertz_mode(Type y, Type mu, Type phi, int give_log) {
  Type one = Type(1.0);
  Type zero = Type(0.0);
  Type log_y = log(y);
  Type log_mu = log(mu);
  Type log_ab = log(phi + one) + phi * log_mu;
  Type log_a = log_ab - log(phi);
  Type log_y_inv_phi = -phi * log_y;
  Type log_y_inv_phi_minus_one = logspace_sub(log_y_inv_phi, zero);
  Type exponential_term = exp(log_a + log_y_inv_phi_minus_one);
  Type log_pdf = log_ab - (phi + one) * log_y - exponential_term;
  
  return give_log ? log_pdf : exp(log_pdf);
}

// ============================================================================
// Quantile: 5 ARGUMENTS (y, mu, phi, tau, give_log)
// ============================================================================

// 9. Kumaraswamy (quantile)
template<class Type>
Type dkum_quantile(Type y, Type mu, Type phi, Type tau, int give_log) {
  Type alpha = log(Type(1.0) - tau) / log(Type(1.0) - pow(mu, phi));
  Type log_res = log(alpha) + log(phi) + (alpha - Type(1.0)) * log(Type(1.0) - pow(y, phi)) + (phi - Type(1.0)) * log(y);
  if(give_log) return log_res; else return exp(log_res);
}

// 10. Vasicek (quantile)
template<class Type>
Type dvasicek_quantile(Type y, Type mu, Type phi, Type tau, int give_log) {
  Type one  = Type(1.0);
  Type half = Type(0.5);
  Type one_minus_phi = one - phi;
  Type z_y   = qnorm(y);
  Type z_mu  = qnorm(mu);
  Type z_tau = qnorm(tau);
  
  Type arg = (sqrt(one_minus_phi) * (z_y - z_mu) + sqrt(phi) * z_tau) / sqrt(phi);
  Type log_pdf = half * (log(one_minus_phi) - log(phi) + z_y * z_y - arg * arg);
  
  return give_log ? log_pdf : exp(log_pdf);
}

// 11. Unit-Weibull (quantile)
template<class Type>
Type dunitweibull(Type y, Type mu, Type phi, Type tau, int give_log)
{
  Type nlog_tau = -log(tau);
  Type nlog_mu  = -log(mu);
  Type nlog_y   = -log(y);
  
  Type log_alpha = log(nlog_tau) - phi * log(nlog_mu);
  Type phi_log_nlog_y = phi * log(nlog_y);
  Type log_pdf = log_alpha + log(phi) - log(y) + 
                 (phi - Type(1.0)) * log(nlog_y) -
                 exp(log_alpha + phi_log_nlog_y);
  
  if(give_log) return log_pdf;
  else return exp(log_pdf);
}

// 12. Unit-Gompertz (quantile)
template<class Type>
Type dunitgompertz_quantile(Type y, Type mu, Type phi, Type tau, int give_log) {
  Type denom = pow(mu, -phi) - Type(1.0);
  Type alpha = (-log(tau) * phi) / denom;
  
  Type term = (alpha / phi) * (pow(y, -phi) - Type(1.0));
  Type log_res = log(alpha) - (phi + Type(1.0)) * log(y) - term;
  if(give_log) return log_res; else return exp(log_res);
}

// 13. Johnson SB (quantile)
template<class Type>
Type djohnsonsb(Type y, Type mu, Type phi, Type tau, int give_log) {
  Type zero = Type(0.0);
  Type one  = Type(1.0);
  
  Type log_y  = log(y);
  Type log_mu = log(mu);
  Type log_one_minus_y =  logspace_sub(zero, log_y);
  Type log_one_minus_mu = logspace_sub(zero, log_mu);
  Type logit_y = log_y - log_one_minus_y;
  Type logit_mu = log_mu - log_one_minus_mu;
  Type alpha = qnorm(tau) - phi * logit_mu;
  Type z = alpha + phi * logit_y;
  Type log_pdf = log(phi) - log_y - log_one_minus_y + dnorm(z, zero, one, 1);
  
  return give_log ? log_pdf : exp(log_pdf);
}

// 14. ASHW (ArcSecant Hyperbolic Weibull) (quantile)
template<class Type>
Type dashw(Type y, Type mu, Type phi, Type tau, int give_log) {
  Type one  = Type(1.0);
  Type half = Type(0.5);
  
  Type h_mu = arcsech(mu);
  Type h_y  = arcsech(y);
  
  Type log_h_mu = log(h_mu);
  Type log_h_y  = log(h_y);
  
  Type log_neg_log_tau = log(-log(tau));
  Type log_alpha = log_neg_log_tau - phi * log_h_mu;
  Type log_term = log_neg_log_tau + phi * (log_h_y - log_h_mu);
  Type term = exp(log_term);
  
  Type log_one_minus_y_squared = log(one - y) + log(one + y);
  Type log_pdf = log_alpha + log(phi) + (phi - one) * log_h_y - term -
                 log(y) - half * log_one_minus_y_squared;
  
  return give_log ? log_pdf : exp(log_pdf);
}

// 15. UBS (Unit-Birnbaum-Saunders) (quantile)
template<class Type>
Type dubs(Type y, Type mu, Type phi, Type tau, int give_log) {
  Type one  = Type(1.0);
  Type two  = Type(2.0);
  Type four = Type(4.0);
  Type half = Type(0.5);
  
  Type log_y   = log(y);
  Type nlog_y  = -log_y;
  Type nlog_mu = -log(mu);
  Type z = -qnorm(tau);
  Type k = phi * z;
  Type root = sqrt(k * k + four);
  Type root_minus_k;
  
  if (asDouble(tau) <= 0.5) {
    root_minus_k = four / (root + k);
  } else {
    root_minus_k = root - k;
  }
  
  Type sqrt_alpha = half * sqrt(nlog_mu) * root_minus_k;
  Type log_alpha = two * log(sqrt_alpha);
  
  Type log_r = half * (log(nlog_y) - log_alpha);
  Type r = exp(log_r);
  Type inverse_r = exp(-log_r);
  Type normal_arg = (r - inverse_r) / phi;
  Type log_term_plus = logspace_add(log_r, -log_r);
  Type log_pdf = dnorm(normal_arg, Type(0.0), one, 1 ) -
                 log_y - log(phi) - log(two) - log(nlog_y) +
                 log_term_plus;
    return give_log ? log_pdf : exp(log_pdf);
}

// ============================================================================
// VECTORIZATION MACROS
// ============================================================================

// Macro for 4 arguments: (y, mu, phi, give_log)
VECTORIZE4_ttti(dbeta_mean)
VECTORIZE4_ttti(dsimplex)
VECTORIZE4_ttti(dvasicek_mean)
VECTORIZE4_ttti(dunitgamma_mean)
VECTORIZE4_ttti(dbessel_mean)
VECTORIZE4_ttti(dbeta_mode)
VECTORIZE4_ttti(dkum_mode)
VECTORIZE4_ttti(dugamma_mode)
VECTORIZE4_ttti(dunitgompertz_mode)
  
// Macro for 5 arguments: (y, mu, phi, tau, give_log)
VECTORIZE5_tttti(dkum_quantile)
VECTORIZE5_tttti(dvasicek_quantile)
VECTORIZE5_tttti(dunitweibull)
VECTORIZE5_tttti(dunitgompertz_quantile)
VECTORIZE5_tttti(djohnsonsb)
VECTORIZE5_tttti(dashw)
VECTORIZE5_tttti(dubs)
  
#endif