library(fungible)
library(noisemaker)
library(future)
library(furrr)

reps <- 100

# Generate population correlation matrix ----------------------------------

mod <- simFA(Model = list(NFac = 1,
                          NItemPerFac = 10),
             Loadings = list(FacLoadDist = "sequential",
                             FacLoadRange = c(.7, .3)))

for (i in 1:reps) {
  Rcb <- noisemaker(mod, method = "CB", target_rmsea = 0.05)
  Rtkl <- noisemaker(mod, method = "TKL", 
                     target_rmsea = Rcb$rmsea,
                     target_cfi = Rcb$cfi,
                     tkl_ctrl = list(penalty = 0))
  
  # cb_fac <- factanal(covmat = Rcb$Sigma, factors = 1, 
  #                    n.obs = 1000, rotation = "none")
  # tkl_fac <- factanal(covmat = Rtkl$Sigma, factors = 1, 
  #                     n.obs = 1000, rotation = "none")
  
  cb_fac_rot <- faMain(R = Rcb$Sigma, n = 1000, numFactors = 1L, 
                       facMethod = "faml", rotate = "targetT",
                       targetMatrix = mod$loadings)$loadings |>  unclass()
  tkl_fac_rot <- faMain(R = Rtkl$Sigma, n = 1000, numFactors = 1L, 
                       facMethod = "faml", rotate = "targetT",
                       targetMatrix = mod$loadings)$loadings |>  unclass()
  
  
}

plan(multisession)