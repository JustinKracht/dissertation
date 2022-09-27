# Check to see whether large minor factor loadings can be "rotated away"
#
# Author: Justin Kracht
# Date: 2022-09-18

library(fungible)
library(noisemaker) # devtools::install_github("JustinKracht/noisemaker")

# Specify a model likely to produce W violations --------------------------

mod <- simFA(Model = list(NFac = 3, NItemPerFac = 5),
             Loadings = list(FacLoadRange = c(.3, .7),
                             FacLoadDist = "sequential"))

# Use the TKL method with a target RMSEA value of 0.08 and no W penalty
# set.seed(42)
tkl_sol <- noisemaker(mod, method = "TKL", target_rmsea = 0.09, 
                  tkl_ctrl = list(penalty = 0))

# How many absolute loadings >= .3 does each minor factor have?
W <- tkl_sol$W # extract the W matrix
apply(W, MARGIN = 2, FUN = function(x) sum(abs(x) >= .3))

# Rotate W to get the largest possible loadings ---------------------------

# Get the principal component loadings of WWt 
X <- tcrossprod(W)
ULU <- eigen(X)
U <- ULU$vectors
L <- diag(sqrt(ULU$values))
P <- U %*% L

out <- svd(W)
Wstar <- varimax(W)$loadings

sol1 <- W %*% t(W)
sol2 <- P %*% t(P)

all.equal(sol1, 
          sol2)

apply(W, MARGIN = 2, FUN = function(x) sum(abs(x) >= .3))
apply(P, MARGIN = 2, FUN = function(x) sum(abs(x) >= .3))

plot(x = 1:15, y = eigen(X)$values)
plot(x = 1:15, y = eigen(X)$values)

# Factor analyze X --------------------------------------------------------



# Rotate the W matrix to a target (null) matrix ---------------------------
# This doesn't converge, but that's ok; it just has to demonstrate that we
# can "rotate away" large minor factor loadings, even if the rotation doesn't
# converge.

target <- matrix(0, nrow = nrow(W), ncol = ncol(W))
rot <- GPArotation::targetQ(L = W, Target = target, maxit = 1e4)

# Check loadings for the first five factors
round(rot$loadings[,1:5], 2)

# How many absolute loadings >= .3 does each minor factor have?
apply(rot$loadings, MARGIN = 2, FUN = function(x) sum(abs(x) >= .3))
