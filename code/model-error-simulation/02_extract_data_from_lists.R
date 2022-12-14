# Extract simulation data from lists and calculate alternative fit statistics

library(dplyr)
library(tidyr)
library(purrr)
library(parallel)
library(stringr)
library(fungible)
library(here)
library(purrrgress)

condition_matrix <- readRDS(here("data", "condition_matrix.RDS"))

# Create results matrix ---------------------------------------------------

reps <- 500

# Create a matrix to hold results
results <- condition_matrix[rep(1:nrow(condition_matrix), each = reps),]

# For each rep in each condition, record the rep number
results$rep_num <- rep(1:reps, times = nrow(condition_matrix))

# Create model error method variable
results <- expand_grid(
  results,
  error_method = c("TKL_rmsea", "TKL_cfi", "TKL_rmsea_cfi", "CB", "WB")
)

# Function definitions ----------------------------------------------------

# Compute the fit statistics that were not computed in the main simulation
# loop

compute_fit_stats <- function(error_method_results, Rpop, k) {
  
  # Check if value is NULL; if so, return empty list
  if (is.null(error_method_results$value)) {
    out = data.frame(RMSEA_thetahat = NA,
                     CFI_thetahat = NA,
                     CRMR_theta = NA,
                     CRMR_thetahat = NA,
                     TLI_theta = NA,
                     TLI_thetahat = NA)
    return(out)
  }
  
  RpopME <- pluck(error_method_results, "value", "Sigma")
  
  p <- ncol(RpopME)
  tp <- p * (p + 1)/2
  
  fout <- factanal(covmat = RpopME, factors = k, 
                   rotation = "none", n.obs = 1000, 
                   control = list(nstart = 100))
  
  Fhat <- fout$loadings
  Sigma_k <- Fhat %*% t(Fhat)
  diag(Sigma_k) <- 1
  
  Tr <- function(X) sum(diag(X))
  
  num_thetahat <- log(det(Sigma_k)) - log(det(RpopME)) + 
    Tr(RpopME %*% solve(Sigma_k)) - p
  
  k <- ncol(Fhat)
  DF <- (p * (p - 1)/2) - (p * k) + (k * (k - 1)/2)
  DF_B <- (p * (p - 1)/2) - (p * p) + (p * (p - 1)/2)
  
  F_T_thetahat <- log(det(Sigma_k)) - log(det(RpopME)) + 
    Tr(RpopME %*% solve(Sigma_k)) - p
  F_B_thetahat <- -log(det(RpopME))
  
  F_T_theta <- log(det(Rpop)) - log(det(RpopME)) + 
    Tr(RpopME %*% solve(Rpop)) - p
  F_B_theta <- -log(det(RpopME))
  
  # Calculate RMSEA thetahat
  RMSEA_thetahat <- sqrt(num_thetahat/DF)
  
  # Calculate CFI thetahat
  CFI_thetahat <- 1 - F_T_thetahat/F_B_thetahat
  
  # Calculate CRMR values
  CRMR_theta <- sqrt(
    sum(((RpopME - Rpop)[upper.tri(Rpop, diag = FALSE)]^2)) / (tp - p)
  )
  
  CRMR_thetahat <- sqrt(
    sum(((RpopME - Sigma_k)[upper.tri(Rpop, diag = FALSE)]^2)) / (tp - p)
  )
  
  # Calculate TLI values (from Xia & Yang, 2019, p. 412)
  TLI_theta <- 1 - ( (F_T_theta / DF) / (F_B_theta / DF_B) )
  TLI_thetahat <- 1 - ( (F_T_thetahat / DF) / (F_B_thetahat / DF_B) )
  
  data.frame(RMSEA_thetahat = RMSEA_thetahat,
             CFI_thetahat = CFI_thetahat,
             CRMR_theta = CRMR_theta,
             CRMR_thetahat = CRMR_thetahat,
             TLI_theta = TLI_theta,
             TLI_thetahat = TLI_thetahat)
}

# Compute d values
compute_d_values <- function(error_method_data, target_rmsea, target_cfi) {
  
  if (is.null(error_method_data$value)) {
    out = data.frame(d1 = NA,
                     d2 = NA,
                     d3 = NA)
    return(out)
  }
  
  rmsea <- pluck(error_method_data, "value", "rmsea")
  cfi <- pluck(error_method_data, "value", "cfi")
  
  d1 <- abs(rmsea - target_rmsea)
  d2 <- abs(cfi - target_cfi)
  d3 <- d1 + d2 
  
  data.frame(
    d1 = d1,
    d2 = d2,
    d3 = d3
  )
}

# Check if there are any major factors in W
check_w_major_factors <- function(W) {
  if (!is.matrix(W)) {
    NA
  } else {
    sum(W[,1] >= .3) > 2
  }
}

# Create a vector of result file paths ------------------------------------

results_files <- list.files(
  path = "data",
  pattern = ".*[0-9]+\\.RDS",
  full.names = TRUE
)

# Read data from each condition and calculate statistics ------------------

# Calculate alternative fit statistics

pbmcapply::pbmclapply(
  X = seq_along(results_files), 
  FUN = function(i, results_files, condition_matrix) {
    
    cat("/nWorking on condition:", i)
    
    # Read in loading matrix
    condition_results <- readRDS(results_files[i])
    
    # Extract the condition number from the results file name
    condition_num <- as.numeric(
      str_extract(results_files[i], pattern = "[0-9]+")
    )
    
    # Extract condition information (loading strength, num. factors, num. items)
    j <- which(condition_matrix$condition_num == condition_num)
    loading_condition <- condition_matrix$loading[j]
    k <- condition_matrix$factors[j]
    p <- condition_matrix$items_per_factor[j] * k
    
    # Extract target RMSEA and CFI values
    target_rmsea <- condition_matrix$target_rmsea[j]
    target_cfi <- condition_matrix$target_cfi[j]
    
    # Create factor loading matrix
    loadings <- switch(loading_condition,
                       "weak" = .4,
                       "moderate" = .6,
                       "strong" = .8)
    
    # Extract Rpop
    mod <- fungible::simFA(
      Model = list(NFac = condition_matrix$factors[j],
                   NItemPerFac = condition_matrix$items_per_factor[j],
                   Model = "oblique"),
      Loadings = list(FacLoadDist = "fixed",
                      FacLoadRange = loadings),
      Phi = list(PhiType = "fixed",
                 MaxAbsPhi = condition_matrix$factor_cor[j])
    )
    
    Rpop <- mod$Rpop
    
    condition_results_matrix <- purrr::map_dfr(
      .x = seq_along(condition_results),
      .f = function(z, Rpop, p, k, target_rmsea, target_cfi,
                    condition_num) {
        
        rep_results <- condition_results[[z]]
        other_fit_stats <- map_dfr(.x = rep_results,
                                   ~ compute_fit_stats(., Rpop = Rpop, k = k))
        d_values <- map_dfr(.x = rep_results,
                            ~ compute_d_values(., 
                                               target_rmsea = target_rmsea, 
                                               target_cfi = target_cfi))
        
        out <- data.frame(
          condition_num = rep(condition_num, 5),
          rep_num = rep(z, 5),
          cfi = map_dbl(rep_results, ~ pluck(., "value", "cfi", .default = NA)),
          cfi_thetahat = other_fit_stats$CFI_thetahat,
          rmsea = map_dbl(rep_results, ~ pluck(., "value", "rmsea", 
                                               .default = NA)),
          rmsea_thetahat = other_fit_stats$RMSEA_thetahat,
          crmr = other_fit_stats$CRMR_theta,
          crmr_thetahat = other_fit_stats$CRMR_thetahat,
          tli = other_fit_stats$TLI_theta,
          tli_thetahat = other_fit_stats$TLI_thetahat,
          m = map_dbl(rep_results, ~ pluck(., "value", "m", .default = NA)),
          v = map_dbl(rep_results, ~ pluck(., "value", "v", .default = NA)),
          eps = map_dbl(rep_results, ~ pluck(., "value", "eps", .default = NA)),
          fn_value = map_dbl(rep_results, ~ pluck(., "value", "fn_value", 
                                                  .default = NA)),
          error_method = c("tkl_rmsea", "tkl_cfi", "tkl_rmsea_cfi", "cb", "wb"),
          w_has_major_factors = map_lgl(rep_results, 
                                        ~ check_w_major_factors(
                                          pluck(., "value", "W", .default = NA)
                                        )),
          d1 = d_values$d1,
          d2 = d_values$d2,
          d3 = d_values$d3,
          warning = map_chr(rep_results, .f = function(x) {
            warning <- pluck(x, "warning")
            if (is.null(warning)) { 
              warning <- NA 
            } else {
              warning <- as.character(warning)
            }
            warning
          }),
          error = map_chr(rep_results, .f = function(x) {
            error <- pluck(x, "error")
            if (is.null(error)) {
              error <- NA
            } else {
              error <- as.character(error)
            }
            error
          })
        )
        
        rownames(out) <- NULL
        out
      },
      Rpop = Rpop, p = p, k = k, target_rmsea = target_rmsea, 
      target_cfi = target_cfi, condition_num
    )
    
    saveRDS(condition_results_matrix, 
            file = paste0("data/results_matrix_", 
                          formatC(i, width = 3, flag = 0), 
                          ".RDS"))
    
  }, results_files = results_files, condition_matrix = condition_matrix, 
  mc.cores = 8
)
