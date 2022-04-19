# Package setup
library(usethis)

# Add a license -----------------------------------------------------------
use_gpl_license(version = 3, include_future = TRUE)

# Add package suggests ----------------------------------------------------
use_package("fungible", type = "Imports")
use_package("MBESS", type = "Imports")
use_package("sem", type = "Imports")
use_package("GA", type = "imports")
use_package("MCMCpack", type = "imports")
use_package("Rdpack", type = "imports")
use_package("MASS", type = "imports")
use_package("crayon", type = "imports")

# Vaccinate gitignore -----------------------------------------------------
git_vaccinate()

# Create a README ---------------------------------------------------------
use_readme_rmd()

# Add testing and test coverage -------------------------------------------
use_testthat()
use_coverage()

# Add badges --------------------------------------------------------------
use_github_actions_badge(name = "test-coverage")
use_lifecycle_badge(stage = "experimental")

# Create a vignette -------------------------------------------------------
use_vignette(
  name = "simulate-model-error",
  title = "Simulating Population Correlation Matrices with Model Error"
)
