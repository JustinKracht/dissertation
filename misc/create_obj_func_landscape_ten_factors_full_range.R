library(furrr)

plan(multisession, workers = 4L)

par_grid <- expand.grid(
  eps = seq(0, 1, length.out = 250),
  nu = seq(0, .75, length.out = 250)
)

simFA_mod <- simFA(Model = list(NFac = 10, 
                                NItemPerFac = 15, 
                                Model = "orthogonal"),
                   Loadings = list(FacLoadRange = .4,
                                   FacLoadDist = "fixed"),
                   ModelError = list(ModelError = TRUE),
                   Seed = 123)

Rpop <- simFA_mod$Rpop
W <- simFA_mod$W

p <- nrow(Rpop)
k <- ncol(simFA_mod$loadings)
df <- (p * (p - 1) / 2) - (p * k) + (k * (k - 1) / 2) # model df
target_rmsea <- 0.09
target_cfi <- 0.9
weights <- list("TKL_RMSEA" = c(1, 0), 
                "TKL_CFI" = c(0, 1), 
                "TKL_RMSEA_CFI" = c(1, 1))
WmaxLoading <- .3
NWmaxLoading <- 2
return_values <- FALSE
penalty <- 1000

CovMajor <- simFA_mod$loadings %*% simFA_mod$Phi %*% t(simFA_mod$loadings)
u <- 1 - diag(CovMajor)

obj_val_list <- vector(mode = "list", length = 3L)
i <- 1
for (error_method_weights in weights) {
  val <- future_map2_dbl(
    .x = par_grid$nu,
    .y = par_grid$eps,
    ~ noisemaker::obj_func(par = c(.x, .y),
                           Rpop = Rpop,
                           W = W,
                           p = p,
                           u = u,
                           df = df,
                           target_rmsea = target_rmsea,
                           target_cfi = target_cfi,
                           weights = weights[[i]],
                           WmaxLoading = WmaxLoading,
                           NWmaxLoading = NWmaxLoading,
                           return_values = return_values,
                           penalty = penalty),
    .progress = TRUE
  )
  obj_val_list[[i]] <- val
  i <- i + 1
}

par_grid <- data.frame(
  eps = rep(par_grid$eps, 3),
  nu  = rep(par_grid$nu, 3),
  error_method = rep(c("TKL[RMSEA]", "TKL[CFI]", "TKL[RMSEA/CFI]"), 
                     each = nrow(par_grid)),
  value = c(obj_val_list[[1]], obj_val_list[[2]], obj_val_list[[3]])
)

par_grid <- par_grid %>%
  mutate(error_method = factor(error_method,
                               levels = c("TKL[RMSEA]", 
                                          "TKL[CFI]", 
                                          "TKL[RMSEA/CFI]")))

without_W_penalty <- par_grid %>%
  filter(value < 1000) %>%
  group_by(error_method) %>%
  mutate(value = (value - min(value, na.rm = TRUE)) / 
           (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) %>%
  ggplot(aes(x = eps, y = nu, fill = value, z = value)) +
  geom_raster() +
  geom_contour(color = "lightgrey") +
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, .15)) +
  facet_grid( ~ error_method, labeller = label_parsed) +
  theme_bw() +
  labs(x = "",
       y = latex2exp::TeX("$\\nu_e$"),
       fill = latex2exp::TeX("$F_{obj}$"),
       title = "Constraints on W Not Violated") +
  theme(panel.spacing = unit(20, "points"))

with_W_penalty <- par_grid %>%
  filter(value >= 1000) %>%
  group_by(error_method) %>%
  mutate(value = (value - min(value, na.rm = TRUE)) / 
           (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) %>%
  ggplot(aes(x = eps, y = nu, fill = value, z = value)) +
  geom_raster() +
  geom_contour(color = "grey") +
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  facet_grid( ~ error_method, labeller = label_parsed) +
  theme_bw() +
  labs(x = "",
       y = latex2exp::TeX("$\\nu_e$"),
       fill = latex2exp::TeX("$F_{obj}$"),
       title = "Constraints on W Violated") +
  theme(panel.spacing = unit(20, "points"))

obj_func_plot <- without_W_penalty / with_W_penalty + 
  plot_layout(guides = "collect")

ggsave(filename = here("img/obj_func_landscape_ten_factor_full_range.pdf"),
       obj_func_plot,
       dpi = 320,
       height = 5.75,
       width = 8)

saveRDS(par_grid, here("data/obj_function_values_ten_factor_full_range.RDS")) 
