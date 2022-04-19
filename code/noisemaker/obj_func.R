#' Objective function for optimizing RMSEA and CFI
#'
#' This is the objective function that is minimized by the `tkl()` function.
#'
#' @param par (vector) Values of model error variance (\eqn{\upsilon}) and
#'   epsilon (\eqn{\epsilon}).
#' @param Rpop (matrix) The model-implied correlation matrix.
#' @param W (matrix) Matrix of provisional minor common factor loadings with
#'   unit column variances.
#' @param p (scalar) Number of variables.
#' @param u (vector) Major common factor variances.
#' @param df (scalar) Model degrees of freedom.
#' @param target_rmsea (scalar) Target RMSEA value.
#' @param target_cfi (scalar) Target CFI value.
#' @param weights (vector) Vector of length two indicating how much weight to
#'   give RMSEA and CFI, e.g., `c(1,1)` (default) gives equal weight to both
#'   indices; `c(1,0)` ignores the CFI value.
#' @param WmaxLoading (scalar) Threshold value for `NWmaxLoading`.
#' @param NWmaxLoading (scalar) Maximum number of absolute loadings \eqn{\ge}
#'   `WmaxLoading` in any column of `W`.
#' @param penalty (scalar) Large (positive) penalty value to apply if the
#'   NWmaxLoading condition is violated.
#' @param return_values (boolean) If `TRUE`, return the objective function value
#'   along with `Rpop`, `RpopME`, `W`, `RMSEA`, `CFI`, `v`, and `eps` values. If
#'   `FALSE`, return only the objective function value.
#'
#' @export

obj_func <- function(par = c(v, eps),
                     Rpop, W, p, u, df,
                     target_rmsea, target_cfi,
                     weights = c(1, 1),
                     WmaxLoading = NULL,
                     NWmaxLoading = 2,
                     penalty = 0,
                     return_values = FALSE) {
  v <- par[1] # error variance
  eps <- par[2] # epsTKL

  # Rescale W using eps
  scaling_matrix <- diag((1 - eps)^(0:(ncol(W) - 1)))
  W <- W %*% scaling_matrix

  # Create W matrix such that the proportion of unique variance accounted for by
  # the minor common factors is v.
  # Adapted from simFA() (lines 691--698)
  wsq <- diag(tcrossprod(W))
  ModelErrorVar <- v * u
  W <- diag(sqrt(ModelErrorVar / wsq)) %*% W
  RpopME <- Rpop + tcrossprod(W)
  diag(RpopME) <- 1

  # ML objective function value for the full model
  # Adapted from simFA() (lines 651--660)
  Ft <- log(det(Rpop)) - log(det(RpopME)) +
    sum(diag(RpopME %*% solve(Rpop))) - p

  # ML objective function value for the baseline (independence) model
  Fb <- -log(det(RpopME))

  # Compute RMSEA and CFI values
  # Adapted from simFA() (lines 651--660)
  rmsea <- sqrt(Ft / df)
  cfi <- 1 - (Ft / -log(det(RpopME)))

  # Define penalty if WmaxLoading and NWmaxLoading are defined
  if (!is.null(WmaxLoading)) {
    # Takes the value 1 if any column of W has more than NWmaxLoading
    # abs(loadings) >= WmaxLoading
    max_loading_indicator <- any(
      max(apply(abs(W) >= WmaxLoading, 2, sum)) > NWmaxLoading
    )
  } else {
    max_loading_indicator <- 0
  }

  weights <- weights / sum(weights) # scale weights to sum to one

  # Compute objective function value
  # fn_value <- weights[1] * (rmsea - target_rmsea)^2 +
  #   weights[2] * (cfi - target_cfi)^2 +
  #   penalty * max_loading_indicator
  fn_value <- weights[1] * ((rmsea - target_rmsea)^2 / target_rmsea^2) +
    weights[2] * (((1 - cfi) - (1 - target_cfi))^2 / (1 - target_cfi)^2) +
    penalty * max_loading_indicator

  # Objective function value weights RMSEA and CFI differences equally; could be
  # changed, if necessary
  if (return_values == FALSE) {
    fn_value
  } else {
    names(fn_value) <- names(v) <- names(eps) <- NULL
    list(
      fn_value = fn_value,
      Rpop = Rpop,
      RpopME = RpopME,
      W = W,
      rmsea = rmsea,
      cfi = cfi,
      v = v,
      eps = eps
    )
  }
}
