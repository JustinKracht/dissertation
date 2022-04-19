#' Generate an sem model from a simFA model object
#'
#' @param mod A `fungible::simFA()` model object.
#'
#' @export
#'
#' @examples
#' ex_mod <- fungible::simFA(Seed = 42)
#' semify(mod = ex_mod)
semify <- function(mod) {
  L <- mod$loadings
  Phi <- mod$Phi
  u <- 1 - mod$h2

  # Specify factor loadings
  loading_spec <- ""
  loadings <- numeric(length = sum(L != 0))
  loading_names <- character(length = sum(L != 0))
  i <- 1
  for (item in seq_len(nrow(L))) {
    for (factor in seq_len(ncol(L))) {
      if (L[item, factor] == 0) next
      loading_spec <- paste0(
        loading_spec, "F", factor, " -> ", "V", item,
        ", lambda", i, ", ", L[item, factor], "\n"
      )
      loadings[i] <- L[item, factor]
      loading_names[i] <- paste0("lambda", i)
      i <- i + 1
    }
  }

  # Specify latent variable variances (fixed at 1)
  latent_var_spec <- ""
  latent_var_names <- colnames(L)
  for (factor in seq_len(ncol(L))) {
    latent_var_spec <- paste0(
      latent_var_spec, "F", factor,
      " <-> ", "F", factor, ", NA, 1\n"
    )
  }

  # Specify latent variable correlations
  latent_cor_spec <- ""
  latent_cor <- NULL
  latent_cor_names <- NULL
  if (ncol(L) > 1) {
    var_pairs <- utils::combn(nrow(Phi), 2)
    latent_cor <- numeric(length = ncol(var_pairs))
    latent_cor_names <- character(length = ncol(var_pairs))
    for (pair in seq_len(ncol(var_pairs))) {
      fi <- var_pairs[1, pair]
      fj <- var_pairs[2, pair]
      latent_cor_spec <- paste0(
        latent_cor_spec, "F", fi, " <-> ", "F", fj,
        ", phi", fi, fj, ", ", Phi[fi, fj], "\n"
      )
      latent_cor[pair] <- Phi[fi, fj]
      latent_cor_names[pair] <- paste0("phi", fi, fj)
    }
  }

  # Specify observed variable variances
  obs_var_spec <- ""
  obs_var <- numeric(length = nrow(L))
  obs_var_names <- paste0("psi", seq_len(nrow(L)))
  for (item in seq_len(nrow(L))) {
    obs_var_spec <- paste0(
      obs_var_spec, "V", item, " <-> ", "V", item, ", psi", item, ", ",
      u[item], "\n"
    )
    obs_var[item] <- u[item]
  }

  # Combine the loading, factor variance, and factor correlation specifications
  # to form the complete model specification
  model <- paste(
    loading_spec,
    latent_var_spec,
    latent_cor_spec,
    obs_var_spec, "\n"
  )

  # Vector of model parameters
  theta <- c(loadings, latent_cor, obs_var)
  names(theta) <- c(loading_names, latent_cor_names, obs_var_names)

  # Return sem model, vector of (named) model parameters, and factor names
  list(
    model = sem::specifyModel(text = model, quiet = TRUE),
    theta = theta,
    latent_var = latent_var_names
  )
}
