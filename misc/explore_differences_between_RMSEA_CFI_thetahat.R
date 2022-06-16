results_matrix %>% 
  filter(rmsea < 3) %>% 
  select(error_method, rmsea, rmsea_thetahat, cfi, cfi_thetahat) %>% 
  mutate(cfi_diff = abs(cfi - cfi_thetahat),
         rmsea_diff = abs(rmsea - rmsea_thetahat)) %>% 
  select(error_method, rmsea_diff, cfi_diff) %>% 
  pivot_longer(-error_method, names_to = "fit_index_diff") %>% 
  group_by(fit_index_diff, error_method) %>% 
  summarise(min = min(value, na.rm = TRUE),
            Q1 = quantile(value, .25, na.rm = TRUE),
            median = quantile(value, .5, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE),
            Q3 = quantile(value, .75, na.rm = TRUE),
            max = max(value, na.rm = TRUE))




results_matrix %>%  
  filter(rmsea < 3) %>% 
  mutate(model_fit_rec = str_replace_all(model_fit_rec, " ", "~"),
         model_fit_rec = fct_relevel(model_fit_rec,
                                     "Fit:~Poor",
                                     "Fit:~Fair",
                                     "Fit:~Very~Good")) %>% 
  ggplot(aes(x = rmsea - rmsea_thetahat)) +
  geom_histogram() +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(limits = c(0, .1)) +
  labs(x = latex2exp::TeX("$RMSEA - RMSEA_{\\hat{\\theta}}$"),
       y = "Count") +
  facet_grid(error_method_rec ~ model_fit_rec,
             labeller = label_parsed) +
  theme_bw(base_size = 9) +
  theme(panel.spacing = unit(8, "pt"))

results_matrix %>%  
  filter(rmsea < 3) %>% 
  mutate(model_fit_rec = str_replace_all(model_fit_rec, " ", "~"),
         model_fit_rec = fct_relevel(model_fit_rec,
                                     "Fit:~Poor",
                                     "Fit:~Fair",
                                     "Fit:~Very~Good")) %>% 
  mutate(error_method_rec = fct_rev(error_method_rec)) %>% 
  ggplot(aes(x = cfi - cfi_thetahat, y = error_method_rec)) +
  geom_boxplot(width = .4, outlier.size = .4, outlier.alpha = .05, size = .4) +
  # scale_x_continuous(limits = c(0, .075)) +
  scale_y_discrete(labels = scales::label_parse()) +
  labs(x = latex2exp::TeX("$CFI - CFI_{\\hat{\\theta}}$"),
       y = "") +
  facet_grid(model_fit_rec * loading_rec ~ factors_rec * items_per_factor_rec,
             labeller = label_parsed) +
  theme_bw()

# EFFECT SIZES

diff_data <- results_matrix %>% 
  mutate(rmsea_diff = rmsea - rmsea_thetahat,
         cfi_diff = cfi_thetahat) %>% 
  mutate(across(all_of(c("factors", "items_per_factor", "loading_numeric", "factor_cor")), ~ scale(.)))

m1 <- lm(rmsea_diff ~ (factors + items_per_factor + loading_numeric + 
                         model_fit + factor_cor + error_method)^2, 
         data = diff_data)


m2 <- lm(cfi_diff ~ (factors + items_per_factor + loading_numeric + 
                         model_fit + factor_cor + error_method)^2, 
         data = diff_data)

coefplot::coefplot(m1, pointSize = 1, color = "black", sort = "magnitude")
coefplot::coefplot(m2, pointSize = 2, color = "black", sort = "magnitude")

# RMSEA DIFF --------------------------------------------------------------

## Factors ----------------------------------------------------------------

results_matrix %>%  
  filter(rmsea < 3) %>% 
  mutate(model_fit_rec = str_replace_all(model_fit_rec, " ", "~"),
         model_fit_rec = fct_relevel(model_fit_rec,
                                     "Fit:~Poor",
                                     "Fit:~Fair",
                                     "Fit:~Very~Good")) %>% 
  mutate(error_method_rec = fct_rev(error_method_rec)) %>% 
  ggplot(aes(x = rmsea - rmsea_thetahat, y = error_method_rec)) +
  geom_boxplot(width = .4, outlier.size = .4, outlier.alpha = .05, size = .4) +
  scale_x_continuous(limits = c(0, .075)) +
  scale_y_discrete(labels = scales::label_parse()) +
  labs(x = latex2exp::TeX("$RMSEA - RMSEA_{\\hat{\\theta}}$"),
       y = "") +
  facet_grid(~ factors_rec,
             labeller = label_parsed) +
  theme_bw()

