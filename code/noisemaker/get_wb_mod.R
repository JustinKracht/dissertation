#' Find an `lm` model to use with the Wu & Browne (2015) model error method
#'
#' The Wu & Browne (2015) model error method takes advantage of the relationship
#' between v and RMSEA:
#'
#' \deqn{v = RMSEA^2 + o(RMSEA^2).}
#'
#' As RMSEA increases, the approximation \eqn{v ~= RMSEA^2} becomes worse. This
#' function generates population correlation matrices with model error for
#' multiple target RMSEA values and then regresses the target RMSEA values on
#' the median observed RMSEA values for each target. The fitted model can then
#' be used to predict a `target_rmsea` value that will give solutions with RMSEA
#' values that are close to the desired value.
#'
#' @param mod A `fungible::simFA()` model object.
#' @param n The number of times to evaluate `wb()` at each point.
#' @param values The number of target RMSEA values to evaluate between 0.02 and
#'   0.1.
#' @param lower (scalar) The smallest target RMSEA value to use.
#' @param upper (scalar) The largest target RMSEA value to use.
#'
#' @return (`lm` object) An `lm` object to use with the \code{\link{wb}}
#'   function to obtain population correlation matrices with model error that
#'   have RMSEA values closer to the target RMSEA values. The `lm` object will
#'   predict a `target_rmsea` value that will give solutions with (median) RMSEA
#'   values close to the desired RMSEA value.
#' @export
#'
#' @examples
#' mod <- fungible::simFA(Seed = 42)
#' set.seed(42)
#' wb_mod <- get_wb_mod(mod)
#' noisemaker(mod, method = "WB", target_rmsea = 0.05, wb_mod = wb_mod)

get_wb_mod <- function(mod, n = 50, values = 10, lower = .01, upper = .095) {
  # Check arguments
  if (!(is.list(mod))) {
    stop("`mod` must be a valid `simFA()` model object.", call. = F)
  }
  if (is.null(mod$loadings) |
      is.null(mod$Phi) |
      is.null(mod$Rpop)) {
    stop("`mod` must be a valid `simFA()` model object.", call. = F)
  }
  if (length(n) != 1L | !is.numeric(n) | n <= 0) {
    stop("`n` must be a number greater than zero.\n",
         crayon::cyan("\u2139"), " You've specified an `n` value of ",
         n, ".", call. = F)
  }
  if (length(values) != 1L | !is.numeric(values) | values < 2) {
    stop("`values` must be a number greater than two.\n",
         crayon::cyan("\u2139"), " You've specified a `values` value of ",
         values, ".", call. = F)
  }
  if (length(lower) != 1L | !is.numeric(lower) | lower <= 0) {
    stop("`lower` must be a number greater than zero.\n",
         crayon::cyan("\u2139"), " You've specified a `lower` value of ",
         lower, ".", call. = F)
  }
  if (length(upper) != 1L | !is.numeric(upper)) {
    stop("`upper` must be a number.\n",
         crayon::cyan("\u2139"), " You've specified an `upper` value of ",
         upper, ".", call. = F)
  }

  k <- ncol(mod$loadings)
  Omega <- mod$Rpop
  p <- nrow(Omega)

  # WB requires m < p; calculate upper bound
  # (1 / target_rmsea^2) > p means that target_rmsea < 1 / sqrt(p)
  max_target_rmsea <- 1 / sqrt(p)
  if (upper >= max_target_rmsea) {
    warning("Specified upper bound was too large and was reduced.")
    # Set new upper bound to the maximum target RMSEA, minus 5% to avoid hitting
    # a boundary
    upper <- max_target_rmsea - 0.05 * max_target_rmsea
  }

  rmsea_values <- seq(lower, upper, length.out = values)

  rmsea_medians <- sapply(
    X = rmsea_values,
    FUN = function(target_rmsea,
                   mod,
                   Omega,
                   k) {
      obs_rmsea <- replicate(n = n, expr = {
        rmsea(wb(mod, target_rmsea, adjust_target = FALSE)$Sigma, Omega, k)
      })
      stats::median(obs_rmsea)
    },
    mod = mod,
    Omega = Omega,
    k = k
  )

  m1 <- stats::lm(rmsea_values ~ poly(rmsea_medians, 2))
  m1
}
