# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

sfpca_arma <- function(X, lambda_u, lambda_v, alpha_u, alpha_v, Omega_u, Omega_v) {
    .Call(`_moma_sfpca_arma`, X, lambda_u, lambda_v, alpha_u, alpha_v, Omega_u, Omega_v)
}

