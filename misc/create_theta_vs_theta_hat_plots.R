rmsea_vs_rmsea_thetahat <- results_matrix %>%  
  filter(rmsea < 3) %>% 
  mutate(model_fit_rec = str_replace_all(model_fit_rec, " ", "~"),
         model_fit_rec = fct_relevel(model_fit_rec,
                                     "Fit:~Poor",
                                     "Fit:~Fair",
                                     "Fit:~Very~Good")) %>% 
  ggplot(aes(x = rmsea, y = rmsea_thetahat)) +
  geom_point(alpha = 0.15, size = .5) +
  geom_abline(slope = 1, intercept = 0)  +
  coord_fixed(xlim = c(0, .3), ylim = c(0, .3)) +
  scale_x_continuous(breaks = c(0, .1, .2, .3)) +
  scale_y_continuous(breaks = c(0, .1, .2, .3)) +
  labs(x = "RMSEA", y = latex2exp::TeX("$RMSEA_{\\hat{\\theta}}$")) +
  facet_grid(model_fit_rec ~ error_method_rec,
             labeller = label_parsed) +
  theme_bw()

ggsave(here("img/rmsea_vs_rmsea_thetahat.png"),
       rmsea_vs_rmsea_thetahat,
       dpi = 320,
       height = 4,
       scale = 1.25)

cfi_vs_cfi_thetahat <- results_matrix %>% 
  mutate(model_fit_rec = str_replace_all(model_fit_rec, " ", "~"),
         model_fit_rec = fct_relevel(model_fit_rec,
                                     "Fit:~Poor",
                                     "Fit:~Fair",
                                     "Fit:~Very~Good")) %>% 
  ggplot(aes(x = cfi, y = cfi_thetahat)) +
  geom_point(alpha = 0.15, size = .5) +
  geom_abline(slope = 1, intercept = 0)  +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, .2, .4, .6, .8, 1)) +
  scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1)) +
  labs(x = "CFI", y = latex2exp::TeX("$CFI_{\\hat{\\theta}}$")) +
  facet_grid(model_fit_rec ~ error_method_rec,
             labeller = label_parsed) +
  theme_bw(base_size = 10)

ggsave(here("img/cfi_vs_cfi_thetahat.png"),
       cfi_vs_cfi_thetahat,
       dpi = 320,
       height = 4,
       scale = 1.25)
