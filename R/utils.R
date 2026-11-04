# utils.R
# Low-level mathematical helpers used throughout the package.

#' Numerically stable sigmoid (logistic) function
#'
#' Computes `1 / (1 + exp(-x))` with input clamping to `[-30, 30]` to prevent
#' overflow or underflow in `exp()`.
#'
#' @param x Numeric vector or matrix of real-valued inputs.
#'
#' @return A numeric object of the same shape as `x` with values in `(0, 1)`.
#'
#' @keywords internal
sigmoid <- function(x) {
  1 / (1 + exp(-pmin(pmax(x, -30), 30)))
}

#' Sample counts from a Negative Binomial distribution (matrix form)
#'
#' Draws integer count matrices from independent Negative Binomial
#' distributions parameterised by mean (`mu`) and dispersion (`theta`).
#' The variance of each entry is `mu + mu^2 / theta`.  Large `theta`
#' approximates a Poisson; small `theta` gives heavy overdispersion.
#'
#' @param mu    Numeric matrix (cells × features) of expected counts.
#'              All values must be non-negative.
#' @param theta Numeric vector of length equal to the number of columns
#'              (features), giving the per-feature dispersion parameter.
#'              Values must be strictly positive.
#'
#' @return An integer matrix of the same dimensions as `mu` containing
#'         NB-distributed counts.
#'
#' @keywords internal
nb_counts <- function(mu, theta) {
  n <- nrow(mu)
  p <- ncol(mu)
  # Broadcast theta (per-feature) to a full n x p matrix
  size_mat <- matrix(theta, nrow = n, ncol = p, byrow = TRUE)
  matrix(
    rnbinom(
      n    = length(mu),
      mu   = as.vector(mu),
      size = as.vector(size_mat)
    ),
    nrow = n,
    ncol = p
  )
}
