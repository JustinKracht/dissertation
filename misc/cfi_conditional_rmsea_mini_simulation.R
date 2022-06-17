library(furrr)

reps <- 100

factors <- unique(results_matrix$factors)
items_per_factor <- unique(results_matrix$items_per_factor)
loading <- unique(results_matrix$loading_numeric)
target_cfi <- seq(0.90, 0.99, by = 0.005)

conditions_matrix <- expand.grid(
  factors = factors,
  items_per_factor = items_per_factor,
  loading = loading,
  target_cfi = target_cfi
)

plan(multisession, workers = 4L)

set.seed(666)
cfi_rmsea_vals <- future_map_dfr(
  .x = 1:nrow(conditions_matrix), 
  .f = function(condition) {
    mod <- simFA(
      Model = list(NFac = conditions_matrix$factors[condition],
                   NItemPerFac = conditions_matrix$items_per_factor[condition],
                   Model = "orthogonal"),
      Loadings = list(FacLoadDist = "fixed",
                      FacLoadRange = conditions_matrix$loading[condition])
    )
    
    map_dfr(.x = 1:reps,
            .f = function(x, target_cfi) {
              sol <- noisemaker(mod = mod, 
                                method = "TKL", 
                                target_rmsea = NULL,
                                target_cfi = target_cfi,
                                tkl_ctrl = list(penalty = 1e6,
                                                NWmaxLoading = 2,
                                                WmaxLoading = .3))
              w_constraints_violated <- sum(sol$W[,1] >= .3) > 2
              c(rmsea = sol$rmsea, cfi = sol$cfi, 
                w_constraints_violated = w_constraints_violated)
            }, target_cfi = conditions_matrix$target_cfi[condition])
  },
  .progress = TRUE
)

cfi_rmsea_vals <- bind_cols(
  cfi_rmsea_vals,
  conditions_matrix[rep(1:nrow(conditions_matrix), each = reps),]
)

cfi_rmsea_vals <- cfi_rmsea_vals %>%
  mutate(
    factors_rec = factor(factors, 
                         labels = levels(results_matrix$factors_rec)),
    loading_rec = factor(loading, 
                         labels = levels(results_matrix$loading_rec)),
    items_per_factor_rec = factor(items_per_factor, 
                                  labels = levels(results_matrix$items_per_factor_rec))
  )

saveRDS(cfi_rmsea_vals, file = here("data/cfi_conditional_rmsea.RDS"))
# rmsea_cfi_vals <- readRDS(here("data/rmsea_conditional_cfi.RDS"))

cfi_rmsea_vals %>%
  ggplot(aes(x = cfi, y = rmsea)) +
  geom_point(size = 1, alpha = .1) +
  geom_abline(aes(slope = -1, intercept = 1), size = .5, alpha = .5) +
  facet_grid(loading_rec * items_per_factor_rec ~ factors_rec) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
  labs(x = "CFI", y = "RMSEA") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(14, units = "points"))

ggsave(filename = here("img/cfi_conditional_rmsea.png"),
       height = 7,
       width = 7,
       scale = 1.15,
       dpi = 320)