# Set the number of reps
reps <- 200

# Create a matrix of fully-crossed conditions
factors <- unique(results_matrix$factors)
items_per_factor <- 15
loading <- c(.4)
target_rmsea <- c(0.090)
penalty <- c(0, .1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)

conditions_matrix <- expand.grid(
  factors = factors,
  items_per_factor = items_per_factor,
  loading = loading,
  target_rmsea = target_rmsea
)


# Function to check whether a W matrix has more than two factor loadings greater
# than 0.3 for any minor factor.
check_constraints <- function(W) {
  any(max(apply(abs(W) >= .3, 2, sum)) > 2)
}

# For each condition, generate a population correlation matrix without model
# error and then generate 200 population correlations with model error using
# the TKL (RMSEA) method
set.seed(666)
constraint_violations <- map_dfr(
  .x = 1:nrow(conditions_matrix), 
  .f = function(condition) {
    
    cat("\nWorking on condition", condition, "of", nrow(conditions_matrix))
    
    # Generate population correlation matrix without model error
    mod <- simFA(
      Model = list(NFac = conditions_matrix$factors[condition],
                   NItemPerFac = conditions_matrix$items_per_factor[condition],
                   Model = "orthogonal"),
      Loadings = list(FacLoadDist = "fixed",
                      FacLoadRange = conditions_matrix$loading[condition])
    )
    
    # Generate 200 population correlation matrices with model error
    pro_map_dfr(.x = 1:reps,
                .f = function(x, target_rmsea) {
                  map_dfr(.x = penalty, 
                          .f = function(penalty, mod, target_rmsea, seed) {
                            set.seed(seed)
                            sol <- noisemaker(mod = mod, 
                                              method = "TKL", 
                                              target_rmsea = target_rmsea,
                                              target_cfi = NULL,
                                              tkl_ctrl = list(penalty = penalty,
                                                              NWmaxLoading = 2,
                                                              WmaxLoading = .3))
                            w_constraints_violated <- check_constraints(sol$W)
                            
                            c(penalty = penalty, 
                              constraints_violated = w_constraints_violated)
                          },
                          mod = mod,
                          target_rmsea = target_rmsea,
                          seed = sample(1e6, 1))
                }, target_rmsea = conditions_matrix$target_rmsea[condition])
  }
)

# Bind the conditions matrix to the results to indicate which result belongs to
# which condition
constraint_violations <- bind_cols(
  conditions_matrix[rep(1:nrow(conditions_matrix), 
                        each = length(penalty) * reps),],
  constraint_violations
)

# Plot the results
constraint_violations %>%
  mutate(factors = as.factor(factors),
         penalty = factor(penalty, 
                          labels = c("0", "0.1", "1", "10", 
                                     "100", "1,000", "10,000",
                                     "100,000", "1,000,000"))) %>%
  group_by(factors, penalty) %>%
  summarise(percent = mean(constraints_violated, na.rm = TRUE)) %>%
  ggplot(aes(x = penalty, y = percent, color = factors, shape = factors, 
             linetype = factors, group = factors)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  scale_color_brewer(palette = "Dark2", type = "qual") +
  labs(y = "Cases with Violtated W Constraints",
       x = latex2exp::TeX("$\\lambda$"),
       color = "Factors", shape = "Factors", linetype = "Factors") +
  theme_bw() +
  theme(legend.position = "bottom")

# Save the plot
ggsave(filename = here("img/penalty_values.png"),
       dpi = 320,
       height = 4,
       width = 6)