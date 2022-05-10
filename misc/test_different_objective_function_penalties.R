
obj_func2 <- function(par = c(v, eps),
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
    cols_that_violate_constraints <- colSums(abs(W) >= .3) > 2
    if (sum(cols_that_violate_constraints) == 0) {
      max_loading_indicator <- 0
    } else {
      max_loading_indicator <- sum(apply(W[cols_that_violate_constraints], 
                                         2, function(w) sum(abs(w[w >= .3]))))
    }
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

tkl2 <- function(mod,
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
            fn = obj_func2,
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
  
  obj_func2(
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

reps <- 20

factors <- unique(results_matrix$factors)
items_per_factor <- 15
loading <- c(.4)
target_rmsea <- c(0.090)
penalty_vec <- c(.1, 1, 10, 100, 1000, 10000)

conditions_matrix <- expand.grid(
  factors = factors,
  items_per_factor = items_per_factor,
  loading = loading,
  target_rmsea = target_rmsea
)

check_constraints <- function(W) {
  any(max(apply(abs(W) >= .3, 2, sum)) > 2)
}

seed_list <- sample(1:1e6, size = nrow(conditions_matrix))

set.seed(666)
constraint_violations <- map_dfr(
  .x = 1:nrow(conditions_matrix), 
  .f = function(condition) {
    
    cat("\nWorking on condition", condition, "of", nrow(conditions_matrix))
    
    mod <- simFA(
      Model = list(NFac = conditions_matrix$factors[condition],
                   NItemPerFac = conditions_matrix$items_per_factor[condition],
                   Model = "orthogonal"),
      Loadings = list(FacLoadDist = "fixed",
                      FacLoadRange = conditions_matrix$loading[condition])
    )
    
    pro_map_dfr(.x = 1:reps,
                .f = function(x, target_rmsea) {
                  map_dfr(.x = penalty_vec, 
                          .f = function(penalty, mod, target_rmsea, seed) {
                            set.seed(seed)
                            sol <- tkl2(mod = mod, 
                                        target_rmsea = target_rmsea,
                                        target_cfi = NULL,
                                        tkl_ctrl = list(penalty = penalty,
                                                        NWmaxLoading = 2,
                                                        WmaxLoading = .3))
                            w_constraints_violated <- check_constraints(sol$W)
                            
                            c(penalty = penalty, 
                              constraints_violated = w_constraints_violated)
                          },
                          mod = mod,
                          target_rmsea = target_rmsea,
                          seed = sample(1e6, 1))
                }, target_rmsea = conditions_matrix$target_rmsea[condition])
  }
)

constraint_violations <- bind_cols(
  conditions_matrix[rep(1:nrow(conditions_matrix), each = length(penalty) * reps),],
  constraint_violations
)

constraint_violations %>%
  mutate(factors = as.factor(factors)) %>%
  group_by(factors, penalty) %>%
  summarise(percent = mean(constraints_violated, na.rm = TRUE),
            penalty = as.factor(penalty)) %>%
  ggplot(aes(x = penalty, y = percent, color = factors, shape = factors, 
             linetype = factors, group = factors)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  scale_color_brewer(palette = "Dark2", type = "qual") +
  labs(y = "Cases with Violtated W Constraints",
       x = latex2exp::TeX("$\\lambda$"),
       color = "Factors", shape = "Factors", linetype = "Factors") +
  theme_bw() +
  theme(legend.position = "bottom")
