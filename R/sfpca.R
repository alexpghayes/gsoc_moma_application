#' @useDynLib moma, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

op_norm <- function(x, A) as.numeric(sqrt(t(x) %*% A %*% x))
norm <- function(x) sqrt(sum(x^2))
prox_l1 <- function(y, lambda) sign(y) * pmax(abs(y) - lambda, 0)

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

  n <- nrow(X)
  p <- ncol(X)

  stopifnot(lambda_u >= 0)
  stopifnot(lambda_v >= 0)
  stopifnot(alpha_u >= 0)
  stopifnot(alpha_v >= 0)
  stopifnot(all(eigen(Omega_u)$values >= 0))  # check if positive semi-definite
  stopifnot(all(eigen(Omega_v)$values >= 0))
  stopifnot(dim(Omega_u) == c(n, n))
  stopifnot(dim(Omega_v) == c(p, p))

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
sfpca_r <- function(X, lambda_u = 0, lambda_v = 0, alpha_u = 0, alpha_v = 0,
                        Omega_u = diag(nrow(X)), Omega_v = diag(ncol(X))) {

  # this function is an R reference implementation

  # some conventions so the following code matches the paper:
  #   - u: left singular vectors, or related to
  #   - v: right singular vectors, or related to
  #   - d: eigenvalues
  #   - vectors are lowercase (u, v), matrices uppercase (X, Omega_u)
  #   - underscores indicate subscripts, LaTeX style
  #   - X is an n x p matrix

  n <- nrow(X)
  p <- ncol(X)

  stopifnot(lambda_u >= 0)
  stopifnot(lambda_v >= 0)
  stopifnot(alpha_u >= 0)
  stopifnot(alpha_v >= 0)
  stopifnot(all(eigen(Omega_u)$values >= 0))  # check if positive semi-definite
  stopifnot(all(eigen(Omega_v)$values >= 0))
  stopifnot(dim(Omega_u) == c(n, n))
  stopifnot(dim(Omega_v) == c(p, p))

  # initialize u, v to rank one SVD solution

  tsvd <- irlba::irlba(X, 1, 1)
  u <- tsvd$u  # n x 1 matrix
  v <- tsvd$v  # p x 1 matrix

  # multiplication by n and p here is to agree with the MATLAB implementation
  S_u <- diag(n) + alpha_u * Omega_u * n
  S_v <- diag(p) + alpha_v * Omega_v * p

  L_u <- max(eigen(S_u)$values) + 0.01
  L_v <- max(eigen(S_v)$values) + 0.01

  tol <- 1e-6
  delta <- 1

  while (delta > tol) {

    delta_u <- delta_v <- 1

    old_u <- u
    old_v <- v

    while (delta_u > tol) {
      tmp_u <- u
      u <- prox_l1(u + (X %*% v - S_u %*% u) / L_u, lambda_u / L_u)
      u <- if (norm(u) > 0) u / op_norm(u, S_u) else u
      delta_u <- norm(u - tmp_u) / norm(tmp_u)
    }

    # same as u case except place with v, and add a single transposition
    while (delta_v > tol) {
      tmp_v <- v
      v <- prox_l1(v + (t(X) %*% u - S_v %*% v) / L_v, lambda_v / L_v)
      v <- if (norm(v) > 0) v / op_norm(v, S_v) else v
      delta_v <- norm(v - tmp_v) / norm(tmp_v)
    }

    delta <- norm(old_u - u) / norm(old_u) + norm(old_v - v) / norm(old_v)
  }

  # normalize singular vectors and calculate singular value
  u <- if (norm(u) > 0) u / norm(u) else u
  v <- if (norm(v) > 0) v / norm(v) else v
  d <- as.numeric(t(u) %*% X %*% v)

  list(u = u, d = d, v = v)
}
