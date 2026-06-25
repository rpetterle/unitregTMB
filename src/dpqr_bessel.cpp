#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// --- Constantes para Integração Gauss-Legendre (20 pontos) ---
const int GL_N = 20;
const double GL_x[] = { -0.9931285991850949, -0.9639719272779138, -0.9122344282513259, -0.8391169718222188, -0.7463319064601508, -0.6360536807265150, -0.5108670019508271, -0.3737060887154195, -0.2277858511416451, -0.0765265211334973, 0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271, 0.6360536807265150, 0.7463319064601508, 0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949 };
const double GL_w[] = { 0.0176140071391521, 0.0406014298003869, 0.0626720483341091, 0.0832767415767048, 0.1019301198172404, 0.1181945319615184, 0.1316886384491766, 0.1420961093183820, 0.1491729864726037, 0.1527533871307258, 0.1527533871307258, 0.1491729864726037, 0.1420961093183820, 0.1316886384491766, 0.1181945319615184, 0.1019301198172404, 0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521 };

// ============================================================================
// PDF, CDF e QF (Escalares)
// ============================================================================

double bessel_pdf_scalar(double y, double mu, double phi) {
  if (y <= 1e-12 || y >= (1.0 - 1e-12)) return 0.0;
  
  double den_zeta = y * (1.0 - y); 
  double diff = y - mu;
  double zeta = std::sqrt(1.0 + (diff * diff) / den_zeta);
  double arg = phi * zeta;
  
  // Bessel K1 ESCALADO
  double k1_scaled = R::bessel_k(arg, 1.0, 2.0);
  if (std::isnan(k1_scaled) || k1_scaled <= 0.0) return 0.0; 
  
  // Cálculo Logarítmico Seguro (Prevenção de Underflow no Denominador)
  double log_correction = phi * (1.0 - zeta);
  double log_num_const = std::log(mu) + std::log(1.0 - mu) + std::log(phi);
  double log_den = std::log(M_PI) + 1.5 * (std::log(y) + std::log(1.0 - y)) + std::log(zeta);
  
  double log_val = log_num_const + std::log(k1_scaled) + log_correction - log_den;
  double result = std::exp(log_val);
  
  if (std::isinf(result) || std::isnan(result)) return 0.0;
  return result;
}

double bessel_cdf_scalar(double q, double mu, double phi) {
  if (q <= 1e-12) return 0.0;
  if (q >= 1.0 - 1e-12) return 1.0;
  
  double half_q = 0.5 * q; 
  double sum = 0.0;
  for (int i = 0; i < GL_N; i++) {
    sum += GL_w[i] * bessel_pdf_scalar(half_q * GL_x[i] + half_q, mu, phi);
  }
  return half_q * sum;
}

double bessel_quantile_scalar(double p, double mu, double phi) {
  if (p <= 0.0) return 0.0;
  if (p >= 1.0) return 1.0;
  
  double lower = 0.0, upper = 1.0, mid = 0.5;
  for(int i = 0; i < 60; ++i) { 
    mid = (lower + upper) * 0.5;
    if(bessel_cdf_scalar(mid, mu, phi) < p) lower = mid; 
    else upper = mid;
  }
  return mid;
}

// ============================================================================
// Gerador Aleatório Inversa Gaussiana Interno (Michael-Schucany-Haas, 1976)
// Remove dependência do pacote statmod
// ============================================================================
double rinvgauss_scalar(double mu, double lambda) {
    double v = R::rnorm(0.0, 1.0);
    double y = v * v;
    double x = mu + (mu * mu * y) / (2.0 * lambda) - 
               (mu / (2.0 * lambda)) * std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y);
    
    double z = R::runif(0.0, 1.0);
    if (z <= (mu / (mu + x))) {
        return x;
    } else {
        return (mu * mu) / x;
    }
}

// ============================================================================
// EXPORTS (Vetorizados)
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_dbessel(NumericVector x, NumericVector mu, NumericVector phi, bool log_prob = false) {
  int n = std::max({(int)x.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double val = bessel_pdf_scalar(x[i % x.size()], mu[i % mu.size()], phi[i % phi.size()]);
    out[i] = log_prob ? (val > 0 ? std::log(val) : R_NegInf) : val;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_pbessel(NumericVector q, NumericVector mu, NumericVector phi, bool lower_tail = true, bool log_prob = false) {
  int n = std::max({(int)q.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double val = bessel_cdf_scalar(q[i % q.size()], mu[i % mu.size()], phi[i % phi.size()]);
    if (!lower_tail) val = 1.0 - val;
    out[i] = log_prob ? std::log(val) : val;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_qbessel(NumericVector p, NumericVector mu, NumericVector phi, bool lower_tail = true, bool log_prob = false) {
  int n = std::max({(int)p.size(), (int)mu.size(), (int)phi.size()});
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double cur_p = p[i % p.size()];
    if (log_prob) cur_p = std::exp(cur_p);
    if (!lower_tail) cur_p = 1.0 - cur_p;
    out[i] = bessel_quantile_scalar(cur_p, mu[i % mu.size()], phi[i % phi.size()]);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rbessel(int n, NumericVector mu, NumericVector phi) {
  NumericVector out(n);
  int nmu = mu.size();
  int nphi = phi.size();
  
  for (int i = 0; i < n; ++i) {
    double cur_mu = mu[i % nmu];
    double cur_phi = phi[i % nphi];
    
    // Parâmetros da decomposição
    double a = cur_mu * cur_phi;
    double b = cur_phi * (1.0 - cur_mu);
    
    // Simula Y1 e Y2 da Inversa Gaussiana
    double y1 = rinvgauss_scalar(a, a * a);
    double y2 = rinvgauss_scalar(b, b * b);
    
    double z = y1 / (y1 + y2);
    
    // Proteção de borda
    if (z <= 1e-12) z = 1e-12;
    if (z >= 1.0 - 1e-12) z = 1.0 - 1e-12;
    
    out[i] = z;
  }
  return out;
}
