#' @useDynLib moma, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

zeros_like <- function(x) rep(0, length(x))  # vectors only
op_norm <- function(x, A) sqrt(t(x) %*% A %*% x)  # operator norm
l2_norm <- function(x) sqrt(sum(x^2))
prox_l1 <- function(y, lambda) sign(y) * pmax(abs(y) - lambda, 0)

max_eigenval <- function(X) {
  # add 0.01 to agree with MATLAB implementation
  max(eigen(X, symmetric = FALSE, only.values = TRUE)$values) + 0.01
}

#' Compute rank 1 sparse and functional principal components
#'
#' @param X A data matrix. Considered n x p by convention.
#' @param lambda_u Sparsity parameter for left singular vectors.
#' @param lambda_v Sparsity parameter for right singular vectors.
#' @param alpha_u Smoothness parameter for left singular vectors.
#' @param alpha_v Smoothness parameter for right singular vectors.
#' @param Omega_u Roughness penalty matrix (n x n) for left singular vectors.
#'   Must be positive semi-definite.
#' @param Omega_v Roughness penalty matrix (p x p) for right singular vectors.
#'   Must be positive semi-definite.
#'
#' @return List with three named elements: left singular vector u, eigenval d,
#'   and right singular vector v.
#' @export
sfpca <- function(X, lambda_u = 0, lambda_v = 0, alpha_u = 0, alpha_v = 0,
                  Omega_u = diag(nrow(X)), Omega_v = diag(ncol(X))) {

  # R wrapper around C++ implementation for lazily argument evaluation and
  # convenient argument checking

  stopifnot(lambda_u >= 0)
  stopifnot(lambda_v >= 0)
  stopifnot(alpha_u >= 0)
  stopifnot(alpha_v >= 0)
  stopifnot(all(svd(Omega_u)$d >= 0))  # check if positive semi-definite
  stopifnot(all(svd(Omega_v)$d >= 0))

  sfpca_arma(X, lambda_u, lambda_v, alpha_u, alpha_v, Omega_u, Omega_v)
}

#' Compute rank 1 sparse and functional principal components
#'
#' @param X A data matrix. Considered n x p by convention.
#' @param lambda_u Sparsity parameter for left singular vectors.
#' @param lambda_v Sparsity parameter for right singular vectors.
#' @param alpha_u Smoothness parameter for left singular vectors.
#' @param alpha_v Smoothness parameter for right singular vectors.
#' @param Omega_u Roughness penalty matrix (n x n) for left singular vectors.
#'   Must be positive semi-definite.
#' @param Omega_v Roughness penalty matrix (p x p) for right singular vectors.
#'   Must be positive semi-definite.
#'
#' @return List with three named elements: left singular vector u, eigenval d,
#'   and right singular vector v.
#' @export
sfpca_r <- function(X, lambda_u = 0, lambda_v = 0, alpha_u = 0, alpha_v = 0,
                        Omega_u = diag(nrow(X)), Omega_v = diag(ncol(X))) {

  stopifnot(lambda_u >= 0)
  stopifnot(lambda_v >= 0)
  stopifnot(alpha_u >= 0)
  stopifnot(alpha_v >= 0)
  stopifnot(all(svd(Omega_u)$d >= 0))  # check if positive semi-definite
  stopifnot(all(svd(Omega_v)$d >= 0))

  # some conventions so the following code matches the paper:
  #   - u: left singular vectors, or related to
  #   - v: right singular vectors, or related to
  #   - d: eigenvalues
  #   - vectors are lowercase (u, v), matrices uppercase (X, Omega_u)
  #   - underscores indicate subscripts, LaTeX style
  #   - X is an n x p matrix

  n <- nrow(X)
  p <- ncol(X)

  # initialize u, v to rank one SVD solution

  tsvd <- irlba::irlba(X, 1, 1)
  u <- tsvd$u  # n x 1 matrix
  v <- tsvd$v  # p x 1 matrix

  # multiplication by n and p here is to agree with the MATLAB implementation
  S_u <- diag(n) + alpha_u * Omega_u * n
  S_v <- diag(p) + alpha_v * Omega_v * p

  L_u <- max_eigenval(S_u)
  L_v <- max_eigenval(S_v)

  tol <- 1e-6
  delta_u <- tol + 1
  delta_v <- tol + 1

  while (delta_u > tol & delta_v > tol) {

    while (delta_u > tol) {
      old_u <- u
      u <- prox_l1(u + (X %*% v - S_u %*% u) / L_u, lambda_u / L_u)
      u_norm <- as.numeric(op_norm(u, S_u))
      u <- if (u_norm > 0) u / u_norm else zeros_like(u)
      delta_u <- l2_norm(u - old_u)
    }

    # same as u case except place with v, and add a single transposition
    while (delta_v > tol) {
      old_v <- v
      v <- prox_l1(v + (t(X) %*% u - S_v %*% v) / L_v, lambda_v / L_v)
      v_norm <- as.numeric(op_norm(v, S_v))
      v <- if (v_norm > 0) v / v_norm else zeros_like(v)
      delta_v <- l2_norm(v - old_v)
    }
  }

  # normalize eigenvectors and calculate eigenvalue
  u <- if (l2_norm(u) > 0) u / l2_norm(u) else u
  v <- if (l2_norm(v) > 0) v / l2_norm(v) else v
  d <- as.numeric(t(u) %*% X %*% v)

  list(u = u, d = d, v = v)
}
