p1 <- rmsea_cfi_vals %>%
  filter(factors %in% c(1, 5, 10),
         loading %in% c(.4, .8),
         items_per_factor == 15) %>% 
  ggplot(aes(x = rmsea, y = cfi)) +
  geom_point(size = 1, alpha = .1) +
  geom_abline(aes(slope = -1, intercept = 1), size = .5, alpha = .5) +
  facet_grid(loading_rec * items_per_factor_rec ~ factors_rec) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
  labs(x = TeX("$RMSEA_{\\Omega}$"), y = TeX("$CFI_{\\Omega}$")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(14, units = "points"))

p2 <- cfi_rmsea_vals %>%
  filter(factors %in% c(1, 5, 10),
         loading %in% c(.4, .8),
         items_per_factor == 15) %>% 
  ggplot(aes(y = rmsea, x = cfi)) +
  geom_point(size = 1, alpha = .1) +
  geom_abline(aes(slope = -1, intercept = 1), size = .5, alpha = .5) +
  facet_grid(loading_rec * items_per_factor_rec ~ factors_rec) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
  labs(y = TeX("$RMSEA_{\\Omega}$"), x = TeX("$CFI_{\\Omega}$")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(14, units = "points"))

library(patchwork)

p3 <- p1 / p2 + plot_annotation(tag_levels = c("A", "B"))

ggsave(filename = "img/conditional_rmsea_and_cfi_plots.png",
       p3,
       dpi = 320,
       height = 6,
       width = 5,
       scale = 1.2)
