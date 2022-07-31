# Cudeck and Browne B matrix
library(fungible)
library(Matrix)
library(assertthat)

# Create symmetric matrix that makes it easy to see what's going on.
S <- matrix(0, nrow = 6, ncol = 6)
S[upper.tri(S, diag = TRUE)] <- 1:(.5 * (6^2 + 6))
S0 <- S
diag(S0) <- 0
S <- S + t(S0)

# Define vec and vecs
vec <- function(A) {
  matrix(as.vector(A), ncol = 1)
}

vecs <- function(A) {
  matrix(as.vector(A)[as.vector(upper.tri(A, diag = TRUE))], ncol = 1)
}

make_transition_matrix <- function(A, sparse = TRUE) {
  p <- nrow(A)
  K <- matrix(0, nrow = p^2, ncol = .5 * (p^2 + p))
  
  # Get the indices of vec(A) that are in the upper triangle of A
  in_upper_tri <- as.vector(upper.tri(A, diag = TRUE))
  in_upper_tri_idx <- which(in_upper_tri)
  
  # Get the indices of vec(A) that are in the lower triangle of A
  in_lower_tri <- as.vector(lower.tri(A, diag = TRUE))
  in_lower_tri_idx <- which(in_lower_tri)
  
  # Remove diagonal elements from the lower triangle so that they're not double-
  # counted
  in_diag <- as.vector(diag(p)[lower.tri(diag(p), diag = TRUE)])
  in_diag_idx <- which(in_diag == 1)
  
  # For each column in T, leave all elements 0 except for the cell corresponding
  # to the index of the next element in vec(A) that is in the upper triangle. 
  # Set that element equal to 1.
  for (i in 1:length(in_upper_tri_idx)) {
    K[in_upper_tri_idx[i],i] <- 1
    K[in_lower_tri_idx[i],i] <- 1
  }
  
  # Off diagonal elements appear twice, so scale them by half
  in_off_diag <- matrix(upper.tri(A, diag = FALSE), ncol = 1) * .5
  
  # The transition matrix will have many zero elements, so using a sparse matrix
  # will speed up computations and conserve memory.
  if (sparse) {
    K <- Matrix::Matrix(K, sparse = TRUE)
  }
  
  return(K)
}

# Skip all the transition matrix BS and create D directly
create_D <- function(p, sparse = TRUE) {
  pstar <- .5 * (10^2 + 10)
  diag_vec <- rep(2, length.out = pstar)
  
  count <- 0
  i <- 0
  while (count < (pstar - 1)) {
    i <- i + 1
    count <- count + i
    diag_vec[count] <- 1
  }
  
  if (sparse) {
    Matrix::Diagonal(x = diag_vec) 
  } else {
    diag(diag_vec)
  }
}

# Test create D:
X <- matrix(0, 10, 10)
X[upper.tri(X, diag = TRUE)] <- create_D(p = 10)
X

# From MBESS:
duplication.matrix <- function(n = 1) {
  if ((n < 1) | (round(n) != n)) 
    stop("n must be a positive integer")
  d <- matrix(0, n * n, n * (n + 1)/2)
  count = 0
  for (j in 1:n) {
    d[(j - 1) * n + j, count + j] = 1
    if (j < n) {
      for (i in (j + 1):n) {
        d[(j - 1) * n + i, count + i] <- 1
        d[(i - 1) * n + j, count + i] <- 1
      }
    }
    count = count + n - j
  }
  return(d)
}

Ktest <- Matrix(duplication.matrix(n = 6), sparse = TRUE)

# Test the transition matrix function
K <- make_transition_matrix(S)
assertthat::are_equal(crossprod(K, vec(S)), vecs(S))

# Create a function that gives the Moore-Penrose Inverse of the K matrix
make_transition_inv <- function(K) {
  solve(crossprod(K)) %*% t(K)
}

Kinv <- make_transition_inv(K)
assertthat::are_equal(crossprod(Kinv, vecs(S)), vec(S))
