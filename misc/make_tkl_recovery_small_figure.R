tkl_recovery <- noisy_data %>%
  filter(factor_corr != "Factor Cor.: 0.3",
         factors %in% c("Factors: 3", "Factors: 10"),
         factor_loading != "Loading: 0.6") %>% 
  ggplot(aes(x = rmsea, y = cfi)) +
  geom_point(alpha = .3, size = 1) +
  geom_hline(aes(yintercept = target_cfi), data = target_values %>%
               filter(factor_corr != "Factor Cor.: 0.3",
                      factors %in% c("Factors: 3", "Factors: 10"),
                      factor_loading != "Loading: 0.6"),
             size = .25) +
  geom_vline(aes(xintercept = target_rmsea), data = target_values %>%
               filter(factor_corr != "Factor Cor.: 0.3",
                      factors %in% c("Factors: 3", "Factors: 10"),
                      factor_loading != "Loading: 0.6"),
             size = .25) +
  scale_x_continuous(n.breaks = 3) +
  facet_grid(factors * factor_loading ~ items_per_factor * factor_corr) +
  labs(x = TeX("$RMSEA_{\\Omega}$"), y = TeX("$CFI_{\\Omega}$")) +
  theme_bw()

ggsave(filename = "img/tkl_recovery_small.png",
       tkl_recovery,
       dpi = 320,
       height = 4.5,
       width = 5,
       scale = 1.15)
