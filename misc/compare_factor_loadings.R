library(fungible)
library(noisemaker)
library(future)
library(furrr)
library(tidyverse)

reps <- 100

# Generate population correlation matrix ----------------------------------

mod <- simFA(Model = list(NFac = 1,
                          NItemPerFac = 10),
             Loadings = list(FacLoadDist = "sequential",
                             FacLoadRange = c(.7, .3)))

out <- purrr::map_dfr(.x = 1:reps, .f = function(i) {
  Rcb <- noisemaker(mod, method = "CB", target_rmsea = 0.08)
  Rtkl <- noisemaker(mod, method = "TKL", 
                     target_rmsea = Rcb$rmsea,
                     target_cfi = Rcb$cfi,
                     tkl_ctrl = list(penalty = 0))
  
  cb_fac_rot <- faMain(R = Rcb$Sigma, n = 1000, numFactors = 1L, 
                       facMethod = "faml", rotate = "targetT",
                       targetMatrix = mod$loadings)$loadings |>  unclass()
  tkl_fac_rot <- faMain(R = Rtkl$Sigma, n = 1000, numFactors = 1L, 
                        facMethod = "faml", rotate = "targetT",
                        targetMatrix = mod$loadings)$loadings |>  unclass()
  
  tibble(
    rep = rep(i, times = 10),
    item = 1:10,
    tkl_rmsea = rep(Rtkl$rmsea, 10),
    tkl_cfi = rep(Rtkl$cfi, 10),
    cb_rmsea = rep(Rcb$rmsea, 10),
    cb_cfi = rep(Rcb$cfi, 10),
    tkl_loadings = tkl_fac_rot[,1],
    cb_loadings = cb_fac_rot[,1]
  )
})

out <- dplyr::mutate(out, loading_diff = tkl_loadings - cb_loadings)

ggplot(out, aes(x = item, y = loading_diff)) +
  geom_jitter() +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  theme_bw()
