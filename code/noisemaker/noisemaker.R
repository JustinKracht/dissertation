#' Simulate a population correlation matrix with model error
#'
#' This tool lets the user generate a population correlation matrix with model
#' error using one of three methods: (1) the Tucker, Koopman, and Linn (TKL;
#' 1969) method, (2) the Cudeck and Browne (CB; 1992) method, or (3) the Wu and
#' Browne (WB; 2015) method. If the CB or WB methods are used, the user can
#' specify the desired RMSEA value. If the TKL method is used, an optimization
#' procedure finds a solution that produces RMSEA and/or CFI values that are
#' close to the user-specified values.
#'
#' @param mod A `fungible::simFA()` model object.
#' @param method (character) Model error method to use ("TKL", "CB", or "WB").
#' @param target_rmsea (scalar) Target RMSEA value.
#' @param target_cfi (scalar) Target CFI value.
#' @param tkl_ctrl (list) A control list containing the following TKL-specific
#'   arguments. See the `tkl()` help file for more details.
#' @param wb_mod (`lm` object) An optional `lm` object used to find a target
#'   RMSEA value that results in solutions with RMSEA values close to the
#'   desired value. Note that if no `wb_mod` is provided, a model will be
#'   estimated at run time. If many population correlation matrices are going to
#'   be simulated using the same model, it will be considerably faster to
#'   estimate `wb_mod` ahead of time. See also `get_wb_mod()`.
#'
#' @return A list containing \eqn{\Sigma}, RMSEA and CFI values, and the TKL
#'   parameters (if applicable).
#' @export
#'
#' @examples
#' mod <- fungible::simFA(Seed = 42)
#'
#' set.seed(42)
#' # Simulate a population correlation matrix using the TKL method with target
#' # RMSEA and CFI values specified.
#' noisemaker(mod, method = "TKL",
#'            target_rmsea = 0.05,
#'            target_cfi = 0.95,
#'            tkl_ctrl = list(optim_type = "optim"))
#'
#' # Simulate a population correlation matrix using the CB method with target
#' # RMSEA value specified.
#' noisemaker(mod, method = "CB",
#'            target_rmsea = 0.05)
#'
#' # Simulation a population correlation matrix using the WB method with target
#' # RMSEA value specified.
#' noisemaker(mod,
#'            method = "WB",
#'            target_rmsea = 0.05)
noisemaker <- function(mod,
                       method = c("TKL", "CB", "WB"),
                       target_rmsea = 0.05,
                       target_cfi = NULL,
                       tkl_ctrl = list(),
                       wb_mod = NULL) {

  if (is.null(target_rmsea) & is.null(target_cfi)) {
    stop("Either target RMSEA or target CFI must be specified.",
         call. = F)
  }
  if (!is.numeric(target_rmsea) & !is.null(target_rmsea)) {
    stop("Target RMSEA value must be a number or NULL.\n",
         crayon::cyan("\u2139"), " You've specified a target RMSEA value of ",
         target_rmsea, ".", call. = F)
  }
  if (!is.numeric(target_cfi) & !is.null(target_cfi)) {
    stop("Target CFI value must be either a number or NULL.\n",
         crayon::cyan("\u2139"), " You've specified a target CFI value of ",
         target_cfi, ".", call. = F)
  }
  if (!is.null(target_rmsea)) {
    if (target_rmsea < 0 | target_rmsea > 1) {
      stop("The target RMSEA value must be a number between 0 and 1.\n",
           crayon::cyan("\u2139"), " You've specified a target RMSEA value of ",
           target_rmsea, ".", call. = F)
    }
  }
  if (!is.null(target_cfi) & (method != "TKL")) {
    stop(
      "The TKL method must be used when a CFI value is specified.\n",
      crayon::cyan("\u2139")," You've selected the ", method," method.\n",
      crayon::cyan("\u2139")," You've specified a target CFI value of ",
      target_cfi, "."
    )
  }
  if (!(method %in% c("TKL", "WB", "CB"))) {
    stop("`method` must be `TKL`, `CB`, or `WB`.\n",
         crayon::cyan("\u2139"), " You've specified ",
         method, " as `method`.", call. = F)
  }
  if (!(is.list(mod)) |
      is.null(mod$loadings) |
      is.null(mod$Phi) |
      is.null(mod$Rpop)) {
    stop("`mod` must be a valid `simFA()` model object.", call. = F)
  }
  if (mod$cn$ModelError$ModelError == TRUE) {
    warning(paste0("The `simFA()` object you provided includes model error",
                   " parameters that will be ignored by this function."))
  }

  out_list <- list(Sigma = NA,
                   rmsea = NA,
                   cfi = NA,
                   fn_value = NA,
                   m = NA,
                   v = NA,
                   eps = NA,
                   W = NA)

  k <- ncol(mod$loadings) # number of major factors

  if (method == "WB") {
    wb_out <- wb(mod = mod,
              target_rmsea = target_rmsea,
              wb_mod = wb_mod)
    out_list$Sigma <- wb_out$Sigma
    out_list$rmsea <- rmsea(out_list$Sigma, mod$Rpop, k)
    out_list$cfi <- cfi(out_list$Sigma, mod$Rpop)
    out_list$m <- wb_out$m
  } else if (method == "CB") {
    out_list$Sigma <- cb(mod = mod,
                         target_rmsea = target_rmsea)
    out_list$rmsea <- rmsea(out_list$Sigma, mod$Rpop, k)
    out_list$cfi <- cfi(out_list$Sigma, mod$Rpop)
  } else if (method == "TKL") {
    tkl_out <- tkl(mod = mod,
                   target_rmsea = target_rmsea,
                   target_cfi = target_cfi,
                   tkl_ctrl = tkl_ctrl)

    out_list$Sigma <- tkl_out$RpopME
    out_list$rmsea <- tkl_out$rmsea
    out_list$cfi <- tkl_out$cfi
    out_list$v <- tkl_out$v
    out_list$eps <- tkl_out$eps
    out_list$W <- tkl_out$W
    out_list$fn_value <- tkl_out$fn_value
  }

  out_list
}
