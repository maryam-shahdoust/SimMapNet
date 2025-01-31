#' Bayesian Estimation of Sparse Precision Matrix
#'
#' This function estimates the precision matrix in a Bayesian framework utilizing additional information about feature similarities and applies sparsity via quantile-based thresholding.
#'
#' @param Y Numeric matrix (n x p): Data matrix with n samples and p features.
#' @param nue Numeric: Prior degrees of freedom for the Wishart distribution.
#' @param distance Numeric matrix (p x p): Distance matrix between features.
#' @param epsilon1 Numeric: Small positive value to ensure positive definiteness in omega.
#' @param epsilon2 Numeric: Small positive value to ensure positive definiteness in final precision matrix.
#' @param alpha Numeric: Kernel width parameter.
#' @param kernel.id Integer (1 or 2): Kernel type (1 for Gaussian, 2 for exponential).
#' @param quantile_level Numeric (0-1): Quantile threshold for sparsification.
#' 
#' @return A list containing:
#' \item{Inverse_sigma_sparsed}{Sparse precision matrix}
#' \item{precision_matrix}{Binary adjacency matrix}
#' 
#' @examples
#' Y <- matrix(rnorm(100), nrow=10, ncol=10)
#' distance <- as.matrix(dist(Y))
#' result <- SimMapNet(Y, nue=5, distance, epsilon1=0.01, epsilon2=0.01, alpha=1, kernel.id=1, quantile_level=0.9)
#' 
#' @export
SimMapNet <- function(Y, nue, distance, epsilon1, epsilon2, alpha, kernel.id, quantile_level) { 
  p <- ncol(Y)  # Number of features
  n <- nrow(Y)  # Number of samples 
  
  # Kernel function
  K <- switch(kernel.id,
              `1` = exp(- (distance) ^ 2 / (2 * alpha^2)),  # Gaussian kernel
              `2` = exp((-1) * distance / alpha),  # Exponential kernel
              stop("Invalid kernel.id. Choose 1 for Gaussian kernel or 2 for Laplacian kernel."))
  
  # Predefine structure for omega according to smartPCA 
  S <- cov(Y)  # Covariance matrix of samples
  V <- sqrt(diag(S))  # Standard deviations
  V_diag <- diag(V)  # Convert to diagonal matrix
  
  # Construct omega.hat
  omega.hat <- V_diag %*% K %*% V_diag
  omega.hat <- omega.hat + diag(epsilon1, p)  # Ensure positive definiteness
  
  # Bayesian estimation for precision matrix
  omega.star <- ((n / (n + nue)) * S) + ((nue / (n + nue)) * omega.hat)
  nue.star <- nue + n
  
  IB <- as.matrix(nue.star * omega.star) + epsilon2  # Ensure positive definiteness
  Isigma_star <- solve(IB)  # Inverse covariance matrix
  
  # Apply sparsification using quantile thresholding
  upper_tri_values <- abs(Isigma_star[upper.tri(Isigma_star)])
  ths <- quantile(upper_tri_values, quantile_level)
  
  if (quantile_level == 0) {
    ths <- ths - 0.0001
  } else if (quantile_level == 1) {
    ths <- ths + 0.0001
  }
  
  # Apply thresholding
  Isigma_sparse <- ifelse(abs(Isigma_star) <= ths, 0, Isigma_star)
  
  # Create binary adjacency matrix
  theta <- ifelse(Isigma_sparse != 0, 1, 0)
  diag(theta) <- 1
  
  return(list(Inverse_sigma_sparsed = Isigma_sparse, precision_matrix = theta))
}
