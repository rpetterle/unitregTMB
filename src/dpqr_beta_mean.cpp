#include <Rcpp.h>
#include <algorithm> // Para std::max

using namespace Rcpp;

// ============================================================================
// Beta Distribution (Mean Parameterization)
// Parameters:
//   mu  : Mean (0 < mu < 1)
//   phi : Precision (phi > 0)
//   shape1 (alpha) = mu * phi
//   shape2 (beta)  = (1 - mu) * phi
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_dbeta_mean(const NumericVector x, 
                             const NumericVector mu, 
                             const NumericVector phi, 
                             const bool log_prob = false) {
    int nx = x.size();
    int nmu = mu.size();
    int nphi = phi.size();
    int n = std::max({nx, nmu, nphi});
    
    NumericVector out(n);

    for (int i = 0; i < n; ++i) {
        double cur_x   = x[i % nx];
        double cur_mu  = mu[i % nmu];
        double cur_phi = phi[i % nphi];

        if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) {
            out[i] = R_NaN;
            continue;
        }

        double shape1 = cur_mu * cur_phi;
        double shape2 = (1.0 - cur_mu) * cur_phi;

        out[i] = R::dbeta(cur_x, shape1, shape2, log_prob ? 1 : 0);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_pbeta_mean(const NumericVector q, 
                             const NumericVector mu, 
                             const NumericVector phi, 
                             const bool lower_tail = true, 
                             const bool log_prob = false) {
    int nq = q.size();
    int nmu = mu.size();
    int nphi = phi.size();
    int n = std::max({nq, nmu, nphi});
    
    NumericVector out(n);

    for (int i = 0; i < n; ++i) {
        double cur_q   = q[i % nq];
        double cur_mu  = mu[i % nmu];
        double cur_phi = phi[i % nphi];

        if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) {
            out[i] = R_NaN;
            continue;
        }

        double shape1 = cur_mu * cur_phi;
        double shape2 = (1.0 - cur_mu) * cur_phi;

        out[i] = R::pbeta(cur_q, shape1, shape2, lower_tail ? 1 : 0, log_prob ? 1 : 0);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_qbeta_mean(const NumericVector p, 
                             const NumericVector mu, 
                             const NumericVector phi, 
                             const bool lower_tail = true, 
                             const bool log_prob = false) {
    int np = p.size();
    int nmu = mu.size();
    int nphi = phi.size();
    int n = std::max({np, nmu, nphi});
    
    NumericVector out(n);

    for (int i = 0; i < n; ++i) {
        double cur_p   = p[i % np];
        double cur_mu  = mu[i % nmu];
        double cur_phi = phi[i % nphi];

        if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) {
            out[i] = R_NaN;
            continue;
        }

        double shape1 = cur_mu * cur_phi;
        double shape2 = (1.0 - cur_mu) * cur_phi;

        out[i] = R::qbeta(cur_p, shape1, shape2, lower_tail ? 1 : 0, log_prob ? 1 : 0);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_rbeta_mean(const int n, 
                             const NumericVector mu, 
                             const NumericVector phi) {
    // Para random generation, o tamanho 'n' é fixo pelo usuário
    NumericVector out(n);
    int nmu = mu.size();
    int nphi = phi.size();

    for (int i = 0; i < n; ++i) {
        double cur_mu  = mu[i % nmu];
        double cur_phi = phi[i % nphi];

        if (cur_mu <= 0.0 || cur_mu >= 1.0 || cur_phi <= 0.0) {
            out[i] = R_NaN;
            continue;
        }

        double shape1 = cur_mu * cur_phi;
        double shape2 = (1.0 - cur_mu) * cur_phi;

        out[i] = R::rbeta(shape1, shape2);
    }
    return out;
}
