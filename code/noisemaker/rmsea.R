#' Calculate RMSEA between two correlation matrices
#'
#' Given two correlation matrices of the same dimension, calculate the RMSEA
#' value using the degrees of freedom for the exploratory factor analysis model
#' (see details).
#'
#' @param Sigma (matrix) Population correlation or covariance matrix (with model
#'   error).
#' @param Omega (matrix) Model-implied population correlation or covariance
#'   matrix.
#' @param k (scalar) Number of major common factors.
#'
#' @details Note that this function uses the degrees of freedom for an
#'   exploratory factor analysis model: \deqn{df = p(p-1)/2-(pk)+k(k-1)/2,}
#'   where \eqn{p} is the number of items and \eqn{k} is the number of major
#'   factors.
#'
#' @md
#' @export
#'
#' @examples
#' mod <- fungible::simFA(Model = list(NFac = 3),
#'                        Seed = 42)
#' set.seed(42)
#' Omega <- mod$Rpop
#' Sigma <- noisemaker(
#'   mod = mod,
#'   method = "CB",
#'   target_rmsea = 0.05
#' )$Sigma
#' rmsea(Sigma, Omega, k = 3)

rmsea <- function(Sigma, Omega, k) {
  if (!is.matrix(Sigma) | !is.matrix(Omega)) {
    stop("Sigma and Omega must be matrices.", call. = F)
  } else if (all.equal(Sigma, t(Sigma)) != TRUE) {
    stop("Sigma must be a symmetric matrix.", call. = F)
  } else if (all.equal(Omega, t(Omega)) != TRUE) {
    stop("Omega must be a symmetric matrix.", call. = F)
  } else if (all.equal(dim(Omega), dim(Sigma)) != TRUE) {
    stop("Sigma and Omega must have the same dimensions.", call. = F)
  } else if (!is.numeric(k) | ((k %% 1) != 0) | k < 0) {
    stop("`k` must be a non-negative integer.\n",
         crayon::cyan("\u2139"), " You've specified ", k, " as `k`.", call. = F)
  }

  p <- nrow(Sigma)
  df <- (p * (p - 1) / 2) - (p * k) + (k * (k - 1) / 2)
  Fm <- log(det(Omega)) - log(det(Sigma)) + sum(diag(Sigma %*% solve(Omega))) - p
  sqrt(Fm / df)
}
