# Show why CFI degrades as the dimensions of Sigma increase
factors <- 1:20
out <- cbind("factors" = factors, "Ft" = NA, "Fb" = NA, "cfi" = NA, "rmsea" = NA)

out <- lapply(
  X = factors, 
  FUN = function(factor_num) {
    m1 <- simFA(Model = list(NFac = factor_num, NItemPerFac = 15),
                Loadings = list(FacLoadRange = .4,
                                FacLoadDist = "fixed"))
    
    error_mod <- noisemaker(mod = m1, method = "TKL", target_rmsea = 0.09)
    
    Omega <- m1$Rpop
    Sigma <- error_mod$Sigma
    p <- nrow(Sigma)
    Ft <- log(det(Omega)) - log(det(Sigma)) + sum(diag(Sigma %*% solve(Omega))) - p
    Fb <- -log(det(Sigma))
    
    c(factors = factor_num, 
      Ft = Ft, 
      Fb = Fb, 
      cfi = cfi(Sigma, Omega), 
      rmsea = rmsea(Sigma, Omega, k = factor_num))
  }
)

out <- out %>% 
  bind_rows() %>% 
  pivot_longer(-factors, names_to = "func_type") %>% 
  mutate(func_type = case_when(func_type == "Ft" ~ "Hypothesized",
                               func_type == "Fb" ~ "Baseline", 
                               TRUE ~ func_type))

p1 <- out %>% 
  filter(func_type %in% c("Hypothesized", "Baseline")) %>% 
  ggplot(aes(x = factors, y = value, color = func_type, linetype = func_type)) + 
  geom_point(size = 1) +
  geom_line() +
  scale_color_brewer(palette = "Dark2", type = "qual") +
  labs(x = "Factors", y = "Minimized Discrepancy Function Value", 
       color = "Model", linetype = "Model") +
  theme_bw() +
  theme(legend.position = "bottom")

p2 <- out %>% 
  filter(func_type == "cfi") %>% 
  ggplot(aes(x = factors, y = value)) + 
  geom_point(size = 1) +
  geom_line() +
  labs(y = "CFI", x = "Factors") +
  theme_bw()

p3 <- out %>% 
  filter(func_type == "rmsea") %>% 
  ggplot(aes(x = factors, y = value)) + 
  geom_point(size = 1) +
  geom_line() +
  scale_y_continuous(limits = c(0, 0.1)) +
  labs(y = "RMSEA", x = "Factors") +
  theme_bw()

cfi_factors_plot <- p1 + (p2 / p3)
