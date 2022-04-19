#' Optimize TKL parameters to find a solution with target RMSEA and CFI values
#'
#' Find the optimal W matrix such that the RMSEA and CFI values are as close as
#' possible to the user-specified target values.
#'
#' @param mod A `fungible::simFA()` model object.
#' @param target_rmsea (scalar) Target RMSEA value.
#' @param target_cfi (scalar) Target CFI value.
#' @param tkl_ctrl (list) A control list containing the following TKL-specific
#'   arguments:
#'   * weights (vector) Vector of length two indicating how much weight to give
#'   RMSEA and CFI, e.g., `c(1,1)` (default) gives equal weight
#'   to both indices; `c(1,0)` ignores the CFI value.
#'   * v_start (scalar) Starting value to use for \eqn{\upsilon}, the proportion
#'   of uniqueness variance reallocated to the minor common factors. Note that
#'   only `v` as a proportion of the unique (not total) variance is supported
#'   in this function.
#'   * eps_start (scalar) Starting value to use for \eqn{\epsilon}, which
#'   controls how common variance is distributed among the minor common factors.
#'   * NminorFac (scalar) Number of minor common factors.
#'   * WmaxLoading (scalar) Threshold value for `NWmaxLoading`.
#'   * NWmaxLoading (scalar) Maximum number of absolute loadings \eqn{\ge}
#'   `WmaxLoading` in any column of \eqn{W}.
#'   * penalty (scalar) Penalty applied to objective function if the
#'   `NmaxLoading` condition isn't satisfied.
#'   * optim_type (character)  Which optimization function to use, `optim` or
#'   `ga`? `optim()` is faster, but might not converge in some cases. If `optim`
#'    doesn't converge, `ga` will be used as a fallback option.
#'   * max_tries (numeric) How many times to restart optimization with new start
#'   parameter values if optimization doesn't converge?
#'   * factr (numeric) controls the convergence of the "L-BFGS-B" method.
#'   Convergence occurs when the reduction in the objective is within this
#'   factor of the machine tolerance. Default is 1e7, that is a tolerance of
#'   about 1e-8. (when using `optim`).
#'   * maxit (number) Maximum number of iterations to use (when using `optim`).
#'   * ncores (boolean/scalar) Controls whether `ga()` optimization is done in
#'   parallel. If `TRUE`, uses the maximum available number of processor cores.
#'   If `FALSE`, does not use parallel processing. If an integer is provided,
#'   that's how many processor cores will be used (if available).
#' @md
#'
#' @details This function attempts to find optimal values of the TKL parameters
#'   \eqn{\upsilon} and \eqn{\epsilon} such that the resulting correlation
#'   matrix with model error (\eqn{\Sigma}) has population RMSEA and/or CFI
#'   values that are close to the user-specified values. It is important to note
#'   that solutions are not guaranteed to produce RMSEA and CFI values that are
#'   reasonably close to the target values; in fact, some combinations of RMSEA
#'   and CFI will be difficult or impossible to obtain for certain models (see
#'   Lai & Green, 2016). It can be particularly difficult to find good solutions
#'   when additional restrictions are placed on the minor factor loadings (i.e.,
#'   using the `WmaxLoading` and `NWmaxLoading` arguments).
#'
#'   Optimization is fastest when the `optim_type = optim` optimization method
#'   is chosen. This indicates that optimization should be done using the
#'   `L-BFGS-B` algorithm implemented in the `optim()` function. However, this
#'   method can sometimes fail to find a solution. In that case, I recommend
#'   setting `optim_type = ga`, which indicates that a genetic algorithm
#'   (implemented in `GA::ga()`) will be used. This method takes longer than
#'   `optim()` but is more likely to find a solution.
#'
#' @export
#' @references Tucker, L. R., Koopman, R. F., & Linn, R. L. (1969). Evaluation
#'   of factor analytic research procedures by means of simulated correlation
#'   matrices. *Psychometrika*, *34*(4), 421–459.
#'   <https://doi.org/10/chcxvf>

tkl <- function(mod,
                target_rmsea = NULL,
                target_cfi = NULL,
                tkl_ctrl = list()) {

  # Create default tkl_ctrl list; modify elements if changed by the user
  tkl_ctrl_default <- list(weights = c(rmsea = 1, cfi = 1),
                           v_start = stats::runif(1, 0.02, 0.9),
                           eps_start = stats::runif(1, 0, 0.8),
                           NMinorFac = 50,
                           WmaxLoading = NULL,
                           NWmaxLoading = 2,
                           debug = FALSE,
                           penalty = 1e6,
                           optim_type = "optim",
                           max_tries = 100,
                           factr = 1e6,
                           maxit = 5000,
                           ncores = FALSE)

  # Update the elements of the default tkl_ctrl list that have been changed by
  # the user
  tkl_ctrl_default <- tkl_ctrl_default[sort(names(tkl_ctrl_default))]
  tkl_ctrl_default[names(tkl_ctrl)] <- tkl_ctrl

  # Create objects for each of the elements in tkl_ctrl
  weights <- tkl_ctrl_default$weights
  v_start <- tkl_ctrl_default$v_start
  eps_start <- tkl_ctrl_default$eps_start
  NMinorFac <- tkl_ctrl_default$NMinorFac
  WmaxLoading <- tkl_ctrl_default$WmaxLoading
  NWmaxLoading <- tkl_ctrl_default$NWmaxLoading
  debug <- tkl_ctrl_default$debug
  penalty <- tkl_ctrl_default$penalty
  optim_type <- tkl_ctrl_default$optim_type
  ncores <- tkl_ctrl_default$ncores
  max_tries <- tkl_ctrl_default$max_tries
  factr <- tkl_ctrl_default$factr

  # Check arguments
  if (!is.null(target_rmsea)) {
    if (target_rmsea < 0 | target_rmsea > 1) {
      stop("The target RMSEA value must be a number between 0 and 1.\n",
           crayon::cyan("\u2139"), " You've specified a target RMSEA value of ",
           target_rmsea, ".", call. = F)
    }
  }
  if (!is.null(target_cfi)) {
    if (target_cfi > 1 | target_cfi < 0) {
      stop("Target CFI value must be between 0 and 1\n",
           crayon::cyan("\u2139"), " You've specified a target CFI value of ",
           target_cfi, ".", call. = F)
    }
  }
  if (is.null(target_cfi) & is.null(target_rmsea)) {
    stop("Either target RMSEA or target CFI (or both) must be specified.")
  }
  if (eps_start < 0 | eps_start > 1) {
    stop("The value of eps_start must be between 0 and 1.", call. = F)
  }
  if (v_start < 0 | v_start > 1) {
    stop("The value of v_start must be between 0 and 1.", call. = F)
  }
  if (!(is.list(mod)) |
      is.null(mod$loadings) |
      is.null(mod$Phi) |
      is.null(mod$Rpop)) {
    stop("`mod` must be a valid `simFA()` model object.", call. = F)
  }
  if (!is.numeric(weights) | length(weights) != 2) {
    stop("`weights` must be a numeric vector of length two.", call. = F)
  }
  if (NMinorFac < 0) {
    stop("The number of minor factors must be non-negative.\n",
         crayon::cyan("\u2139"), " You've asked for ", NMinorFac,
         " minor factors.", call. = F)
  }
  if (!(optim_type %in% c("optim", "ga"))) {
    stop("`optim_type` must be either `optim` or `ga`.\n",
         crayon::cyan("\u2139"), " You've supplied ", optim_type,
         " as `optim_type`.", call. = F)
  }
  if (!is.numeric(penalty) | penalty < 0) {
    stop("`penalty` must be a postive number.\n",
         crayon::cyan("\u2139"), " You've supplied ", penalty,
         " as `penalty`.", call. = F)
  }
  if (!is.null(WmaxLoading)) {
    if (!is.numeric(WmaxLoading) | WmaxLoading <= 0) {
      stop("`WmaxLoading` must be a positive number.\n",
           crayon::cyan("\u2139"), " You've supplied ", WmaxLoading,
           " as `WmaxLoading`.", call. = F)
    }
  }
  if (((NWmaxLoading %% 1) != 0) | NWmaxLoading < 0) {
    stop("`NWmaxLoading` must be a non-negative integer.\n",
         crayon::cyan("\u2139"), " You've supplied ", NWmaxLoading,
         " as `NWmaxLoading`.", call. = F)
  }

  # If no CFI value is given, set the weight to zero and set target_cfi to a
  # no-null value (it will be ignored in the optimization)
  if (is.null(target_cfi)) {
    weights[2] <- 0
    target_cfi <- 999
  }

  # Same for RMSEA
  if (is.null(target_rmsea)) {
    weights[1] <- 0
    target_rmsea <- 999
  }

  L <- mod$loadings
  Phi <- mod$Phi

  # Create W with eps = 0
  W <- MASS::mvrnorm(
    n = nrow(L),
    mu = rep(0, NMinorFac),
    Sigma = diag(NMinorFac)
  )

  p <- nrow(L) # number of items
  k <- ncol(L) # number of major factors

  CovMajor <- L %*% Phi %*% t(L)
  u <- 1 - diag(CovMajor)
  Rpop <- CovMajor
  diag(Rpop) <- 1 # ensure unit diagonal

  df <- (p * (p - 1) / 2) - (p * k) + (k * (k - 1) / 2) # model df

  start_vals <- c(v_start, eps_start)

  if (optim_type == "optim") {
    ctrl <- list(factr = factr)
    if (debug == TRUE) {
      ctrl$trace <- 5
      ctrl$REPORT <- 1
    }
    # Try optim(); if it fails, then use GA instead
    opt <- NULL
    tries <- 0
    converged <- FALSE
    while (converged == FALSE & (tries <= max_tries)) {
      if (tries > 1) start_vals <- c(v_start = stats::runif(1, 0.02, 0.9),
                                     eps_start = stats::runif(1, 0, 0.8))
      tryCatch(
        {
          opt <- stats::optim(
            par = start_vals,
            fn = obj_func,
            method = "L-BFGS-B",
            lower = c(0.001, 0), # can't go lower than zero
            upper = c(1, 1), # can't go higher than one
            Rpop = Rpop,
            W = W,
            p = p,
            u = u,
            df = df,
            target_rmsea = target_rmsea,
            target_cfi = target_cfi,
            weights = weights,
            WmaxLoading = WmaxLoading,
            NWmaxLoading = NWmaxLoading,
            control = ctrl,
            penalty = penalty
          )
          par <- opt$par
        },
        error = function(e) NULL
      )

      tries <- tries + 1
      if (is.null(opt)) {
        converged <- FALSE
        next
      }

      converged <- opt$convergence == 0
    }

    # If the algorithm fails to converge or produces NULL output, try GA instead
    if (is.null(opt)) {
      opt <- list(convergence = FALSE)
    }

    if (opt$convergence != 0) {
      optim_type <- "ga"
      warning("`optim()` failed to converge, using `ga()` instead.",
              call. = FALSE)
    }
  }

  if (optim_type == "ga") {
    opt <- GA::ga(
      type = "real-valued",
      fitness = function(x) {
        -obj_func(x,
          Rpop = Rpop,
          W = W,
          p = p,
          u = u,
          df = df,
          target_rmsea = target_rmsea,
          target_cfi = target_cfi,
          weights = weights,
          WmaxLoading = WmaxLoading,
          NWmaxLoading = NWmaxLoading,
          penalty = penalty
        )
      },
      lower = c(0, 0),
      upper = c(1, 1),
      popSize = 50,
      maxiter = 1000,
      run = 100,
      parallel = ncores,
      monitor = FALSE
    )
    par <- opt@solution[1, ]
  }

  obj_func(
    par = par,
    Rpop = Rpop,
    W = W,
    p = p,
    u = u,
    df = df,
    target_rmsea = target_rmsea,
    target_cfi = target_cfi,
    weights = weights,
    WmaxLoading = WmaxLoading,
    NWmaxLoading = NWmaxLoading,
    return_values = TRUE,
    penalty = penalty
  )
}
