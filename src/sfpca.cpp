// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec prox_l1(arma::vec y, double lambda) {
    // there has to be a better way to do this
    arma::vec z = arma::zeros<arma::vec>(y.n_elem);
    arma::mat both = arma::join_horiz(arma::abs(y) - lambda, z);
    return sign(y) % arma::max(both, 1);
}

double op_norm(arma::vec & x, arma::mat & A) {
    arma::mat quad_form = x.t() * A * x;
    return std::sqrt(quad_form(0, 0));
}

// [[Rcpp::export]]
Rcpp::List sfpca_arma(arma::mat & X, double lambda_u, double lambda_v,
                      double alpha_u, double alpha_v,
                      arma::mat Omega_u, arma::mat Omega_v) {

    // following same notation conventions as in sfpca_r

    int n = X.n_rows;
    int p = X.n_cols;

    // initialize u, v to rank 1 SVD solution

    arma::mat U, V;
    arma::vec s;

    svd(U, s, V, X);

    arma::vec u = U.col(0); // first column of U
    arma::vec v = V.col(0);

    arma::mat I_n, I_p;
    I_n.eye(n, n);
    I_p.eye(p, p);

    // again multiply by n and p to agree with MATLAB implementation
    arma::mat S_u = I_n + alpha_u * Omega_u * n;
    arma::mat S_v = I_p + alpha_v * Omega_v * p;

    // again add fudge factor to the largest eigenval to match MATLAB code
    double L_u = std::real(max(eig_gen(S_u))) + 0.01;
    double L_v = std::real(max(eig_gen(S_v))) + 0.01;

    double tol = 1e-6;
    double delta_u = tol + 1;
    double delta_v = tol + 1;

    while (delta_u > tol && delta_v > tol) {

        while (delta_u > tol) {
            arma::vec old_u = u;
            u = prox_l1(u + (X * v - S_u * u) / L_u, lambda_u / L_u);
            double u_norm = op_norm(u, S_u);
            if (u_norm > 0) {
                u = u / u_norm;
            }
            delta_u = norm(u - old_u);
        }

        while (delta_v > tol) {
          arma::vec old_v = v;
          v = prox_l1(v + (X.t() * u - S_v * v) / L_v, lambda_v / L_v);
          double v_norm = op_norm(v, S_v);
          if (v_norm > 0) {
              v = v / v_norm;
          }
          delta_v = norm(v - old_v);
        }
    }

    double u_norm = norm(u);
    double v_norm = norm(v);

    if (u_norm > 0) {
        u = u / u_norm;
    }

    if (v_norm > 0) {
        v = v / v_norm;
    }

    // can this be done without an intermediate eig matrix?
    arma::mat eig = u.t() * X * v;
    double d = eig(0, 0);

    return Rcpp::List::create(Rcpp::Named("u") = u,
                              Rcpp::Named("d") = d,
                              Rcpp::Named("v") = v);
}

