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

generate_sp_se_cor <- function(mean_sens, mean_spec, 
                               sd_sens, sd_spec, 
                               rho = 0, 
                               n = 1000, is_psa,
                               seed = NULL) {
  #' Generate Correlated Sensitivity and Specificity Samples
  #'
  #' This function generates probabilistic sensitivity and specificity samples
  #' from a correlated bivariate normal distribution on the logit scale.
  #' The Beta-distributed sensitivity and specificity are approximated by
  #' corresponding logit-normal distributions that preserve their means,
  #' variances, and specified correlation.
  #'
  #' @param mean_sens Numeric. Mean sensitivity (on probability scale, 0–1).
  #' @param mean_spec Numeric. Mean specificity (on probability scale, 0–1).
  #' @param sd_sens Numeric. Standard deviation of sensitivity (probability scale).
  #' @param sd_spec Numeric. Standard deviation of specificity (probability scale).
  #' @param rho Numeric. Correlation coefficient between logit(sensitivity) and logit(specificity).
  #' @param n Integer. Number of samples to generate.
  #' @param is_psa Logical. Whether to generate PSA samples.
  #' @param seed Integer (optional). Random seed for reproducibility.
  #' @return A data.frame with columns `sensitivity` and `specificity`.
  #' @examples
  #' generate_sp_se_cor(mean_sens = 0.85, mean_spec = 0.90,
  #'                    sd_sens = 0.05, sd_spec = 0.04,
  #'                    rho = -0.3, n = 1000)
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for this function. Please install it.")
  }
  
  # Optional reproducibility
  if (!is.null(seed)) set.seed(seed)
  
  # Input validation
  if (any(c(mean_sens, mean_spec) <= 0 | c(mean_sens, mean_spec) >= 1))
    stop("Means must be between 0 and 1.")
  if (any(c(sd_sens, sd_spec) <= 0))
    stop("Standard deviations must be positive.")
  if (abs(rho) > 1)
    stop("Correlation (rho) must be between -1 and 1.")
  
  if (is_psa) {
    # Helpers
    logit    <- function(p) log(p / (1 - p))
    inv_logit <- function(x) 1 / (1 + exp(-x))
    
    # Convert mean/sd to alpha, beta parameters for Beta distribution
    calc_ab <- function(mean, sd) {
      var <- sd^2
      tmp <- (mean * (1 - mean)) / var - 1
      alpha <- mean * tmp
      beta  <- (1 - mean) * tmp
      return(list(alpha = alpha, beta = beta))
    }
    
    ab_sens <- calc_ab(mean_sens, sd_sens)
    ab_spec <- calc_ab(mean_spec, sd_spec)
    
    # Exact logit mean and variance using digamma/trigamma
    mean_logit_sens <- digamma(ab_sens$alpha) - digamma(ab_sens$beta)
    mean_logit_spec <- digamma(ab_spec$alpha) - digamma(ab_spec$beta)
    var_logit_sens  <- trigamma(ab_sens$alpha) + trigamma(ab_sens$beta)
    var_logit_spec  <- trigamma(ab_spec$alpha) + trigamma(ab_spec$beta)
    
    # Covariance matrix on logit scale
    sigma <- matrix(c(var_logit_sens, rho * sqrt(var_logit_sens * var_logit_spec),
                      rho * sqrt(var_logit_sens * var_logit_spec), var_logit_spec), 
                    ncol = 2, byrow = TRUE)
    
    # Sample from multivariate normal on logit scale
    mean_logit_vec <- c(mean_logit_sens, mean_logit_spec)
    samples <- mvtnorm::rmvnorm(n = n, mean = mean_logit_vec, sigma = sigma, method = "chol")
    
    # Transform back to probability scale
    df <- data.frame(
      sensitivity = inv_logit(samples[, 1]),
      specificity = inv_logit(samples[, 2])
    )
  } else {
    df <- data.frame(
      sensitivity = mean_sens,
      specificity = mean_spec)
  }
  return(df)
}
