#### f_input helper functions ----

# check function using: docstring(generate_static)
# check function using: docstring(generate_beta)
# check function using: docstring(generate_gamma)
# check function using: docstring(generate_lognormal)
# check function using: docstring(generate_dirichlet)

generate_static <- function(value, n, is_psa) {
  #' Generate Static Values
  #'
  #' This function generates a static value or repeats the value for probabilistic sensitivity analysis (PSA).
  #'
  #' @param value Numeric. The static value to be used.
  #' @param n Integer. The number of simulations (for PSA).
  #' @param is_psa Logical. Whether to generate multiple values for PSA.
  #' @return A numeric vector or a single value depending on `is_psa`.
  #' @examples
  #' generate_static(0.035, 10, TRUE)
  if (is_psa) {
    rep(value, n)
  } else {
    value
  }
}

generate_beta <- function(mean, sigma, n, is_psa, alpha = NULL, beta = NULL) {
  #' Generate Values from a Beta Distribution
  #'
  #' This function generates values from a beta distribution based on mean and standard deviation or
  #' directly provided alpha and beta parameters.
  #'
  #' @param mean Numeric. Mean of the distribution.
  #' @param sigma Numeric. Standard deviation of the distribution.
  #' @param n Integer. Number of samples to generate.
  #' @param is_psa Logical. Whether to generate PSA samples.
  #' @param alpha Numeric (optional). Alpha parameter for the beta distribution.
  #' @param beta Numeric (optional). Beta parameter for the beta distribution.
  #' @return A numeric vector of sampled values or a single mean value.
  #' @examples
  #' generate_beta(0.5, 0.1, 100, TRUE)
  if (!is.null(alpha) & !is.null(beta)) {
    if (is_psa) {
      rbeta(n, alpha, beta)
    } else {
      alpha / (alpha + beta)
    }
  } else {
    variance <- sigma^2
    alpha <- mean * ((mean * (1 - mean)) / variance - 1)
    beta <- (1 - mean) * ((mean * (1 - mean)) / variance - 1)
    
    if (is_psa) {
      rbeta(n, alpha, beta)
    } else {
      mean
    }
  }
}

generate_gamma <- function(mu, sigma, n, is_psa) {
  #' Generate Values from a Gamma Distribution
  #'
  #' This function generates values from a gamma distribution based on mean and standard deviation.
  #'
  #' @param mu Numeric. Mean of the distribution.
  #' @param sigma Numeric. Standard deviation of the distribution.
  #' @param n Integer. Number of samples to generate.
  #' @param is_psa Logical. Whether to generate PSA samples.
  #' @return A numeric vector of sampled values or a single mean value.
  #' @examples
  #' generate_gamma(10, 2, 100, TRUE)
  shape <- (mu / sigma)^2
  scale <- sigma^2 / mu
  
  if (is_psa) {
    rgamma(n, shape = shape, scale = scale)
  } else {
    mu
  }
}

generate_lognormal <- function(mean, ci_low, ci_high, n, is_psa) {
  #' Generate Values from a Log-Normal Distribution
  #'
  #' This function generates values from a log-normal distribution based on mean and confidence interval.
  #'
  #' @param mean Numeric. Mean of the distribution.
  #' @param ci_low Numeric. Lower bound of the confidence interval.
  #' @param ci_high Numeric. Upper bound of the confidence interval.
  #' @param n Integer. Number of samples to generate.
  #' @param is_psa Logical. Whether to generate PSA samples.
  #' @return A numeric vector of sampled values or a single mean value.
  #' @examples
  #' generate_lognormal(1, 0.5, 1.5, 100, TRUE)
  sdlog <- (log(ci_high) - log(ci_low)) / (2 * qnorm(0.975))
  
  if (is_psa) {
    rlnorm(n, meanlog = log(mean), sdlog = sdlog)
  } else {
    mean
  }
}

generate_dirichlet <- function(alpha, n, name, is_psa) {
  #' Generate Samples from a Dirichlet Distribution
  #'
  #' This function generates samples from a Dirichlet distribution for a given alpha vector.
  #'
  #' @param alpha Numeric vector. Parameters of the Dirichlet distribution.
  #' @param n Integer. Number of samples to generate.
  #' @param name Character. Prefix for naming columns of the output data frame.
  #' @param is_psa Logical. Whether to generate PSA samples.
  #' @return A data frame containing sampled values, with columns named using `name` and indices.
  #' @examples
  #' generate_dirichlet(c(1, 2, 3), 10, "prob", TRUE)
  if (is_psa) {
    samples <- matrix(rgamma(n * length(alpha), shape = alpha), nrow = n, byrow = TRUE)
    samples <- samples / rowSums(samples) # Normalize rows to sum to 1
  } else {
    samples <- matrix(alpha / sum(alpha), nrow = 1)
  }
  
  samples_df <- as.data.frame(samples)
  colnames(samples_df) <- paste0(name, "_", seq_along(alpha))
  
  samples_df
}
