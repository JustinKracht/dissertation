# Check to see if RMSEA/CFI combinations can be recovered by noisemaker()

library(fungible)
library(noisemaker)
library(tidyverse)

# Define conditions ----

factors <- c(1, 3, 5, 10)
items_per_factor <- c(5, 15)
factor_corr <- c(0, .3, .6)
factor_loading <- c(.4, .6, .8)

conditions_matrix <- expand.grid(
  factors = factors,
  items_per_factor = items_per_factor,
  factor_corr = factor_corr,
  factor_loading = factor_loading
)

conditions_matrix <- filter(conditions_matrix,
                            !(factors == 1 & factor_corr != 0))

# Generate a simFA model for each condition -------------------------------

simFA_mod_list <- purrrgress::pro_pmap(
  .l = conditions_matrix,
  ~ simFA(Model = list(NFac = ..1,
                       NItemPerFac = ..2,
                       Model = "oblique"),
          Loadings = list(FacLoadRange = ..4),
          Phi = list(MaxAbsPhi = ..3,
                     PhiType = "fixed"),
          ModelError = list(ModelError = TRUE),
          Seed = 123)
)

# Extract RMSEA and CFI fit statistics ----

fit_stats <- purrr::map(
  .x = simFA_mod_list,
  ~ list(RMSEA = pluck(.x, "ModelErrorFitStats", "RMSEA_theta"),
         CFI = pluck(.x, "ModelErrorFitStats", "CFI_theta"))
)

# Use noisemaker to generate models with error ----

set.seed(123)

noisy_mod_list <- purrrgress::pro_pmap(
  .l = list(mod = simFA_mod_list, 
            target_rmsea = purrr::map(fit_stats, ~ pluck(., "RMSEA")),
            target_cfi = purrr::map(fit_stats, ~ pluck(., "CFI"))),
  function(mod, target_rmsea, target_cfi) {
    suppressWarnings(
      replicate(n = 50, 
                noisemaker(
                  mod = mod, method = "TKL", target_rmsea = target_rmsea,
                  target_cfi = target_cfi, 
                  tkl_ctrl = list(factr = 1e7,
                                  max_tries = 1000,
                                  WmaxLoading = .3,
                                  NWmaxLoading = 2,
                                  factr = 1e5))
      )
    )
  }
)

# Unpack list data into a data frame for plotting ----

noisy_data <- map_dfr(noisy_mod_list, 
                      function(x) {as.data.frame(t(x)[,-c(1,5,8)])})
noisy_data <- noisy_data %>%
  mutate(across(rmsea:eps, ~unlist(.)))

noisy_data$condition <- rep(1:nrow(conditions_matrix), each = 50)
noisy_data$target_rmsea <- rep(purrr::map_dbl(fit_stats, ~ pluck(., "RMSEA")),
                               each = 50)
noisy_data$target_cfi <- rep(purrr::map_dbl(fit_stats, ~ pluck(., "CFI")),
                             each = 50)
noisy_data$condition <- as.factor(noisy_data$condition)
noisy_data <- as_tibble(cbind(noisy_data, 
                              conditions_matrix[noisy_data$condition,]))

# Create variable labels
noisy_data <- noisy_data %>%
  mutate(factors = factor(factors,
                          levels = c(1, 3, 5, 10),
                          labels = paste("Factors:", c(1, 3, 5, 10))),
         items_per_factor = factor(items_per_factor,
                                   levels = c(5, 15),
                                   labels = paste("Items/Factor:", c(5, 15))),
         factor_corr = factor(factor_corr,
                              levels = c(0, .3, .6),
                              labels = paste("Factor Cor.:", c(0.0, 0.3, 0.6))),
         factor_loading = factor(factor_loading,
                                 levels = c(0.4, 0.6, 0.8),
                                 labels = paste("Loading:", c(0.4, 0.6, 0.8))))

saveRDS(noisy_data, file = "data/noisy_data.RDS")
