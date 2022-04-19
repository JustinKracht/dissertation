#' Calculate CFI for two correlation matrices
#'
#' Given two correlation matrices of the same dimension, calculate the CFI value
#' value using the independence model as the null model.
#'
#' @param Sigma (matrix) Population correlation or covariance matrix (with model
#'   error).
#' @param Omega (matrix) Model-implied population correlation or covariance
#'   matrix.
#'
#' @export
#'
#' @examples
#' library(fungible)
#' library(noisemaker)
#'
#' mod <- fungible::simFA(Model = list(NFac = 3),
#'                        Seed = 42)
#' set.seed(42)
#' Omega <- mod$Rpop
#' Sigma <- noisemaker(
#'   mod = mod,
#'   method = "CB",
#'   target_rmsea = 0.05
#' )$Sigma
#' cfi(Sigma, Omega)
cfi <- function(Sigma, Omega) {
  if (!is.matrix(Sigma) | !is.matrix(Omega)) {
    stop("Sigma and Omega must be matrices.", call. = F)
  } else if (all.equal(Sigma, t(Sigma)) != TRUE) {
    stop("Sigma must be a symmetric matrix.", call. = F)
  } else if (all.equal(Omega, t(Omega)) != TRUE) {
    stop("Omega must be a symmetric matrix.", call. = F)
  } else if (all.equal(dim(Omega), dim(Sigma)) != TRUE) {
    stop("Sigma and Omega must have the same dimensions.", call. = F)
  }

  p <- nrow(Sigma)
  Ft <- log(det(Omega)) - log(det(Sigma)) +
    sum(diag(Sigma %*% solve(Omega))) - p
  cfi <- 1 - (Ft / -log(det(Sigma)))
  cfi
}
