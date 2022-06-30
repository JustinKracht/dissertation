reps <- 50

factors <- unique(results_matrix$factors)
items_per_factor <- unique(results_matrix$items_per_factor)
loading <- unique(results_matrix$loading_numeric)
factor_cor <- c(0, 0.6)
target_rmsea <- seq(0.025, 0.1, by = 0.005)

conditions_matrix <- expand.grid(
  factors = factors,
  items_per_factor = items_per_factor,
  loading = loading,
  factor_cor = factor_cor,
  target_rmsea = target_rmsea
)

set.seed(666)
rmsea_cfi_vals <- pro_map_dfr(
  .x = 1:nrow(conditions_matrix), 
  .f = function(condition) {
    mod <- simFA(
      Model = list(NFac = conditions_matrix$factors[condition],
                   NItemPerFac = conditions_matrix$items_per_factor[condition],
                   Model = "oblique"),
      Phi = list(MaxAbsPhi = conditions_matrix$factor_cor[condition],
                 PhiType = "fixed"),
      Loadings = list(FacLoadDist = "fixed",
                      FacLoadRange = conditions_matrix$loading[condition])
    )
    
    map_dfr(.x = 1:reps,
            .f = function(x, target_rmsea) {
              sol <- noisemaker(mod = mod, 
                                method = "TKL", 
                                target_rmsea = target_rmsea,
                                target_cfi = NULL,
                                tkl_ctrl = list(penalty = 1e6,
                                                NWmaxLoading = 2,
                                                WmaxLoading = .3))
              w_constraints_violated <- sum(sol$W[,1] >= .3) > 2
              c(rmsea = sol$rmsea, cfi = sol$cfi, 
                w_constraints_violated = w_constraints_violated)
            }, target_rmsea = conditions_matrix$target_rmsea[condition])
  }
)

rmsea_cfi_vals <- bind_cols(
  rmsea_cfi_vals,
  conditions_matrix[rep(1:nrow(conditions_matrix), each = reps),]
)

rmsea_cfi_vals <- rmsea_cfi_vals %>%
  mutate(
    factors_rec = factor(factors, 
                         labels = levels(results_matrix$factors_rec)),
    loading_rec = factor(loading, 
                         labels = levels(results_matrix$loading_rec)),
    items_per_factor_rec = factor(items_per_factor, 
                                  labels = levels(results_matrix$items_per_factor_rec)),
    factor_cor_rec = factor(factor_cor,
                            labels = levels(results_matrix$factor_cor_rec)[c(1,3)])
  )

rmsea_cfi_vals %>%
  ggplot(aes(x = rmsea, y = cfi)) +
  geom_point(size = 1, alpha = .1) +
  geom_abline(aes(slope = -1, intercept = 1), size = .5, alpha = .5) +
  facet_grid(factors_rec * items_per_factor_rec ~ loading_rec * factor_cor_rec) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
  labs(x = "RMSEA", y = "CFI") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(14, units = "points"))

ggsave(here("img/rmsea_conditional_cfi_with_factor_cor.png"),
       height = 10,
       width = 8,
       dpi = 320)
saveRDS(rmsea_cfi_vals, file = here("data/rmsea_cfi_vals_with_factor_cor.RDS"))

# ---------------------------------


x <- map(.x = 1:15,
         .f = function(i) {
           m1 <- simFA(Model = list(NFac = i, NItemPerFac = 5),
                       Loadings = list(FacLoadRange = .4,
                                       FacLoadDist = "fixed"))
           
           set.seed(42)
           error_mod <- noisemaker(mod = m1, method = "TKL", target_rmsea = 0.09)
           
           Sigma <- error_mod$Sigma
           Omega <- m1$Rpop
           diff <- sum(diag(Sigma %*% solve(Omega)))
           
           c("detSigma" = det(Sigma),
             "detOmega" = det(Omega),
             "diff" = diff)
         })

x <- bind_rows(x)
x$factors <- 1:15

x %>% 
  mutate(logDetSigma = log(detSigma),
         logDetOmega = log(detOmega),
         p = factors * 5) %>% 
  select(factors, logDetSigma, logDetOmega, diff, p) %>% 
  pivot_longer(logDetSigma:p, names_to = "type", values_to = "value") %>% 
  ggplot(aes(x = factors, y = value, color = type)) +
  geom_point() +
  geom_line() +
  theme_bw()


f <- 3
pf <- 10
m1 <- simFA(Model = list(NFac = f, NItemPerFac = pf),
            Loadings = list(FacLoadRange = .4,
                            FacLoadDist = "fixed"))
set.seed(123)
error_mod <- noisemaker(mod = m1, method = "TKL", target_rmsea = 0.09)
Sigma <- error_mod$Sigma

X <- m1$Rpop
det(X)
eigen(X)$values

det(Sigma)




df <- function(p, k) {
  (p * (p - 1)/2) - (p * k) + (k * (k - 1)/2)
}

df1 <- df(p = 10, k = 2)
df2 <- df(p = 15, k = 2)

df2 / df1