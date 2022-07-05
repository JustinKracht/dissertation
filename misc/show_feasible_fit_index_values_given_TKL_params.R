library(fungible)

m1 <- simFA(Model = list(NFac = 10,
                         NItemPerFac = 10),
            Loadings = list(FacLoadRange = 0.4,
                            FacLoadDist = "fixed"))

set.seed(666)
sigma_tkl <- noisemaker(m1, method = "TKL", target_rmsea = 0.09)
sigma_wb <- noisemaker(m1, method = "WB", target_rmsea = 0.09)

par(mfrow = c(1,2))
plot(1:nrow(sigma_tkl$Sigma), eigen(sigma_tkl$Sigma)$values, type = "l")
plot(1:nrow(sigma_wb$Sigma), eigen(sigma_wb$Sigma)$values, type = "l")


# Tool to show RMSEA and CFI values from a range of v values --------------

v_range <- c(0.02, .15)
mod <- simFA() # hypothesized model

ex <- map(1:1000,
    .f = function(x, v_range) {
      mod <- simFA(ModelError = list(ModelError = TRUE,
                                     ModelErrorType = "U",
                                     ModelErrorVar = runif(1, min = v_range[1], max = v_range[2]),
                                     epsTKL = runif(1, 0, 1),
                                     WmaxLoading = 1))
    }, v_range = v_range)

fit_stats <- map_dfr(ex, .f = function(x) c("RMSEA" = pluck(x, "ModelErrorFitStats", "RMSEA_thetahat"),
                                        "CFI" = pluck(x, "ModelErrorFitStats", "CFI_thetahat")))
plot(fit_stats$RMSEA, fit_stats$CFI)