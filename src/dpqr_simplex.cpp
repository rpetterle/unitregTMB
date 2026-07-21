#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// Weights and Nodes for Gauss-Legendre (20 points) 
const double GL_X[] = {
    0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271, 0.6360536807265150,
    0.7463319064601508, 0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949
};
const double GL_W[] = {
    0.1527533871307258, 0.1491729864726037, 0.1420961093183820, 0.1316886384491733, 0.1181945319615184,
    0.1019301198172404, 0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521
};

inline double log_pdf_simplex_scalar(double x, double mu, double phi) {
    if (x <= 1e-9 || x >= (1.0 - 1e-9)) return R_NegInf;
    
    double diff = x - mu;
    double one_minus_x = 1.0 - x;
    double one_minus_mu = 1.0 - mu;
    
    // Simplex deviance component
    double d = (diff * diff) / 
               (x * one_minus_x * mu * mu * one_minus_mu * one_minus_mu);
               
    // phi represents sigma^2, the dispersion parameter
    double log_res = -0.5 * (std::log(2.0 * M_PI) + std::log(phi) + 
                     3.0 * (std::log(x) + std::log(one_minus_x))) - 
                     d / (2.0 * phi);
                     
    return log_res;
}

inline double pdf_simplex_scalar(double x, double mu, double phi) {
    return std::exp(log_pdf_simplex_scalar(x, mu, phi));
}

inline double cdf_simplex_scalar(double q, double mu, double phi) {
    if (q <= 0.0) return 0.0;
    if (q >= 1.0) return 1.0;

    double a = 0.0;
    double b = q;
    double center = 0.5 * (b + a);
    double half_width = 0.5 * (b - a);
    double sum = 0.0;

    for (int i = 0; i < 10; ++i) {
        double dx = half_width * GL_X[i];
        double val_plus = pdf_simplex_scalar(center + dx, mu, phi);
        double val_minus = pdf_simplex_scalar(center - dx, mu, phi);
        sum += GL_W[i] * (val_plus + val_minus);
    }

    double cdf = half_width * sum;
    if (cdf < 0.0) return 0.0;
    if (cdf > 1.0) return 1.0;
    return cdf;
}

inline double invcdf_simplex_scalar(double p, double mu, double phi) {
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return 1.0;
    
    double lower = 0.0, upper = 1.0;
    double mid = 0.5;
    
    for(int i = 0; i < 60; ++i) {
        mid = 0.5 * (lower + upper);
        double val = cdf_simplex_scalar(mid, mu, phi);
        if (val < p) lower = mid;
        else upper = mid;
    }
    return mid;
}

// [[Rcpp::export]]
NumericVector cpp_dsimplex(const NumericVector x, const NumericVector mu, const NumericVector phi, const bool log_prob = false) {
    int n = std::max({(int)x.size(), (int)mu.size(), (int)phi.size()});
    NumericVector out(n);
    int nx = x.size(), nmu = mu.size(), nphi = phi.size();

    for (int i = 0; i < n; ++i) {
        double cur_mu = mu[i % nmu];
        double cur_phi = phi[i % nphi];
        
        if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0) {
            out[i] = R_NaN;
        } else {
            double cur_x = x[i % nx];
            if (cur_x <= 0 || cur_x >= 1) {
                out[i] = log_prob ? R_NegInf : 0.0;
            } else {
                double val = log_pdf_simplex_scalar(cur_x, cur_mu, cur_phi);
                out[i] = log_prob ? val : std::exp(val);
            }
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_psimplex(const NumericVector q, const NumericVector mu, const NumericVector phi, const bool lower_tail = true, const bool log_prob = false) {
    int n = std::max({(int)q.size(), (int)mu.size(), (int)phi.size()});
    NumericVector out(n);
    int nq = q.size(), nmu = mu.size(), nphi = phi.size();

    for (int i = 0; i < n; ++i) {
        double cur_mu = mu[i % nmu];
        double cur_phi = phi[i % nphi];
        
        if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0) {
            out[i] = R_NaN;
        } else {
            double val = cdf_simplex_scalar(q[i % nq], cur_mu, cur_phi);
            if (!lower_tail) val = 1.0 - val;
            out[i] = log_prob ? std::log(val) : val;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_qsimplex(const NumericVector p, const NumericVector mu, const NumericVector phi, const bool lower_tail = true, const bool log_prob = false) {
    int n = std::max({(int)p.size(), (int)mu.size(), (int)phi.size()});
    NumericVector out(n);
    int np = p.size(), nmu = mu.size(), nphi = phi.size();

    for (int i = 0; i < n; ++i) {
        double cur_mu = mu[i % nmu];
        double cur_phi = phi[i % nphi];
        
        if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0) {
            out[i] = R_NaN;
            continue;
        }

        double cur_p = p[i % np];
        if (log_prob) cur_p = std::exp(cur_p);
        if (!lower_tail) cur_p = 1.0 - cur_p;
        
        if (cur_p < 0 || cur_p > 1) out[i] = R_NaN;
        else out[i] = invcdf_simplex_scalar(cur_p, cur_mu, cur_phi);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_rsimplex(const int n, const NumericVector mu, const NumericVector phi) {
    NumericVector out(n);
    NumericVector u = runif(n);
    int nmu = mu.size(), nphi = phi.size();

    for (int i = 0; i < n; ++i) {
        double cur_mu = mu[i % nmu];
        double cur_phi = phi[i % nphi];
        
        if (cur_mu <= 0 || cur_mu >= 1 || cur_phi <= 0) {
            out[i] = R_NaN;
        } else {
            out[i] = invcdf_simplex_scalar(u[i], cur_mu, cur_phi);
        }
    }
    return out;
}
