#' Wu & Browne model error method
#'
#' Generate a population correlation matrix using the model described in Wu and
#' Browne (2015).
#'
#' @param mod A `fungible::simFA()` model object.
#' @param target_rmsea (scalar) Target RMSEA value.
#' @param wb_mod (`lm` object) An optional `lm` object used to find a target
#'   RMSEA value that results in solutions with RMSEA values close to the
#'   desired value. Note that if no `wb_mod` is provided, a model will be
#'   estimated at run time. If many population correlation matrices are going to
#'   be simulated using the same model, it will be considerably faster to
#'   estimate `wb_mod` ahead of time. See also `get_wb_mod()`.
#' @param adjust_target (TRUE; logical) Should the target_rmsea value be
#'   adjusted to ensure that solutions have RMSEA values that are close to the
#'   provided target RMSEA value? Defaults to TRUE and should stay there unless
#'   you have a compelling reason to change it.
#'
#' @author Justin Kracht <krach018@umn.edu>
#' @references Wu, H., & Browne, M. W. (2015). Quantifying adventitious error in
#'   a covariance structure as a random effect. *Psychometrika*, *80*(3),
#'   571–600. <https://doi.org/10/gjrkc4>
#'
#' @export
#' @details The Wu and Browne method generates a correlation matrix with model
#'   error (\eqn{\Sigma}) using
#'
#'   \deqn{(\Sigma | \Omega) ~ IW(m, m \Omega),}
#'
#'   where \eqn{m ~= 1/\epsilon^2} is a precision parameter related to RMSEA
#'   (\eqn{\epsilon}) and \eqn{IW(m, m \Omega)} denotes an inverse Wishart
#'   distribution. Note that *there is no guarantee that the RMSEA will be very
#'   close to the target RMSEA*, particularly when the target RMSEA value is
#'   large. Based on experience, the method tends to give solutions with RMSEA
#'   values that are larger than the target RMSEA values. Therefore, it might be
#'   worth using a target RMSEA value that is somewhat lower than what is
#'   actually needed. Alternatively, the \code{\link{get_wb_mod}} function can
#'   be used to estimate a coefficient to shrink the target RMSEA value by an
#'   appropriate amount so that the solution RMSEA values are close to the
#'   (nominal) target values.
#'
#' @examples
#' # Specify a default model using simFA()
#' mod <- fungible::simFA(Seed = 42)
#'
#' set.seed(42)
#' wb(mod, target_rmsea = 0.05)

wb <- function(mod,
               target_rmsea,
               wb_mod = NULL,
               adjust_target = TRUE) {

  if (!(is.list(mod)) |
      is.null(mod$loadings) |
      is.null(mod$Phi) |
      is.null(mod$Rpop)) {
    stop("`mod` must be a valid `simFA()` model object.", call. = F)
  }
  if (target_rmsea < 0 | target_rmsea > 1) {
    stop("The target RMSEA value must be a number between 0 and 1.\n",
         crayon::cyan("\u2139"), " You've specified a target RMSEA value of ",
         target_rmsea, ".", call. = F)
  }
  if (!is.null(wb_mod)) {
    if (class(wb_mod) != "lm") {
      stop("`wb_mod` must be an object of class `lm`.", call. = F)
    }
  }

  if (is.null(wb_mod) & (adjust_target == TRUE)) {
    if (target_rmsea >= 0.095) {
      upper <- target_rmsea + 0.01
    } else {
      upper <- 0.095
    }
    wb_mod <- get_wb_mod(mod, upper = upper)
  }

  # Use wb_mod to find the correct target_rmsea value to use
  if (!is.null(wb_mod) & (adjust_target == TRUE)) {
    target_rmsea <- stats::predict(
      wb_mod,
      newdata = data.frame(rmsea_medians = target_rmsea)
    )
  }

  v <- target_rmsea^2
  m <- v^-1 # m is the precision parameter, Wu and Browne (2015), p. 576

  Omega <- mod$Rpop
  p <- nrow(Omega)
  if (m < p) {
    stop("Target RMSEA value is too large, try a smaller value.", call. = FALSE)
  }

  Sigma <- MCMCpack::riwish(m, m * Omega)
  list(Sigma = stats::cov2cor(Sigma),
       m = m)
}
