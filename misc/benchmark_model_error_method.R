library(microbenchmark)

mod <- simFA(Model = list(NFac = 10, NItemPerFac = 5), 
             Loadings = list(FacLoadRange = .4, FacLoadDist = "fixed"),
             Seed = 123)
safe_noisemaker <- possibly(.f = noisemaker, otherwise = NA, quiet = TRUE)

times <- microbenchmark(
  safe_noisemaker(mod, method = "TKL", target_rmsea = 0.05),
  safe_noisemaker(mod, method = "TKL", target_rmsea = 0.05, target_cfi = 0.95),
  safe_noisemaker(mod, method = "TKL", target_rmsea = NULL, target_cfi = 0.95),
  safe_noisemaker(mod, method = "WB", target_rmsea = 0.05),
  safe_noisemaker(mod, method = "CB", target_rmsea = 0.05),
  times = 10
)
