log_sum_exp <- function(log_vals) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::logSumExp(log_vals)
  } else {
    m <- max(log_vals)
    m + log(sum(exp(log_vals - m)))
  }
}

# No underrun condition
NURS_stop <- function(log_vals, log_eps_h) {
  max(log_vals[1], log_vals[length(log_vals)]) <=
    log_eps_h + log_sum_exp(log_vals)
}

# Recursive sub stopping
NURS_sub_stop <- function(log_vals, log_eps_h) {
  n <- length(log_vals)
  if (n < 2) return(FALSE)
  if (NURS_stop(log_vals, log_eps_h)) return(TRUE)
  left_indices <- 1:(n / 2)
  NURS_sub_stop(log_vals[left_indices], log_eps_h) ||
    NURS_sub_stop(log_vals[-left_indices], log_eps_h)
}

NURS_step <- function(logpdf, theta, epsilon, h, M) {
  d <- length(theta)
  # hit
  z <- rnorm(d)
  rho <- z / sqrt(sum(z^2))

  # run
  s <- runif(1, -h / 2, h / 2)
  log_cur <- logpdf(theta)
  log_shifted <- logpdf(theta + s * rho)
  theta0 <- if (
    (log_shifted - log_cur >= 0) || (runif(1) < exp(log_shifted - log_cur))
  )
    theta + s * rho else theta

  # helpers
  log_eps_h <- log(epsilon) + log(h)
  orbit_points <- list(theta0)
  log_vals <- logpdf(theta0)

  # bookkeeping to get ends
  left <- right <- 1

  # Orbit selection procedure
  B <- sample(c(FALSE, TRUE), M, replace = TRUE)
  for (k in seq_len(M)) {
    B_k <- B[k]
    n_ext <- 2^(k - 1)
    orbit_ext <- lapply(
      seq_len(n_ext),
      \(i)
        (if (B_k) orbit_points[[right]] else orbit_points[[left]]) +
          (2 * B_k - 1) * i * h * rho
    )
    log_ext <- sapply(orbit_ext, logpdf)

    if (NURS_sub_stop(log_ext, log_eps_h)) break

    if (B_k) {
      orbit_points <- c(orbit_points, orbit_ext)
      log_vals <- c(log_vals, log_ext)
    } else {
      orbit_points <- c(orbit_ext, orbit_points)
      log_vals <- c(log_ext, log_vals)
      left <- 1
    }
    right <- right + n_ext

    if (NURS_stop(log_vals, log_eps_h)) break
  }

  orbit_points[[sample(
    length(log_vals),
    1,
    prob = exp(log_vals - log_sum_exp(log_vals))
  )]]
}

NURS <- function(logpdf, theta_init, n, epsilon, h, M) {
  d <- length(theta_init)
  draws <- matrix(NA_real_, n, d)
  draws[1, ] <- theta_init
  for (i in 2:n) {
    draws[i, ] <- NURS_step(logpdf, draws[i - 1, ], epsilon, h, M)
  }
  draws
}
