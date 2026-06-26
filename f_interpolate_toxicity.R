f_interpolate_toxicity <- function(prevalence_data , interpolation_method = "linear", 
                                   monthly_cycles = 0:60, plot_output = FALSE) {
  #' Interpolate toxicity prevalence over time
  #'
  #' Interpolates observed toxicity prevalence measured at discrete follow-up
  #' time points to estimate monthly toxicity prevalence over the model time
  #' horizon.
  #'
  #' The interpolated prevalence curves are used as inputs to the state-
  #' transition model to estimate toxicity-related costs, disutilities and the
  #' impact of preventative interventions over time.
  #'
  #' @param prevalence_data  Data frame containing observed toxicity prevalence
  #'   with the following columns:
  #'   \itemize{
  #'     \item `month` — Follow-up time in months.
  #'     \item `toxicity` — Observed prevalence of the toxicity.
  #'   }
  #' @param interpolation_method Character string specifying the interpolation
  #'   method. Supported options are `"linear"` (default) and `"spline"`.
  #' @param monthly_cycles Numeric vector specifying the time points (months) at
  #'   which interpolated prevalence should be estimated. Defaults to `0:60`.
  #' @param plot_output Logical. If `TRUE`, plots the observed and interpolated
  #'   toxicity prevalence curves. Default is `FALSE`.
  #'
  #' @details
  #' Clinical studies typically report toxicity prevalence at a limited number
  #' of follow-up visits (e.g., baseline, 2, 12, 24, 36 and 60 months). This
  #' function interpolates these observations to obtain monthly prevalence
  #' estimates required by the cohort state-transition model.
  #'
  #' Linear interpolation is used by default because it provides a transparent
  #' approximation between observed data points. Spline interpolation can be
  #' used for exploratory analyses or model validation.
  #'
  #' @return A data frame with two columns:
  #' \itemize{
  #'   \item `month` — Model cycle (months).
  #'   \item `toxicity` — Interpolated toxicity prevalence.
  #' }
  #'
  #' @examples
  #' prevalence_data  <- data.frame(
  #'   month = c(0, 2, 12, 24, 36, 60),
  #'   toxicity = c(0.05, 0.08, 0.12, 0.18, 0.20, 0.22)
  #' )
  #'
  #' toxicity_curve <- f_interpolate_toxicity(
  #'   prevalence_data  = prevalence_data ,
  #'   interpolation_method = "linear",
  #'   plot_output = TRUE
  #' )
  #'
  #' head(toxicity_curve)
  #'
  #' @export
  
  # Helper function: Validate inputs
  validate_inputs <- function(prevalence_data , interpolation_method) {
    if (!is.data.frame(prevalence_data ) || !all(c("month", "toxicity") %in% names(prevalence_data ))) {
      stop("`prevalence_data ` must be a data frame with columns `month` and `toxicity`.")
    }
    if (!interpolation_method %in% c("linear", "spline")) {
      stop("`interpolation_method` must be either 'linear' or 'spline'.")
    }
  }
  
  # Helper function: Perform interpolation
  perform_interpolation <- function(method, prevalence_data , monthly_cycles) {
    if (method == "linear") {
      interpolated <- approx(
        x = prevalence_data $month, 
        y = prevalence_data $toxicity, 
        xout = monthly_cycles
      )
    } else if (method == "spline") {
      interpolated <- spline(
        x = prevalence_data $month, 
        y = prevalence_data $toxicity, 
        xout = monthly_cycles
      )
    }
    data.frame(month = interpolated$x, toxicity = interpolated$y)
  }
  
  # Helper function: Plot results
  plot_results <- function(original_data, interpolated_data, method) {
    plot(
      original_data$month, original_data$toxicity, 
      type = "p", pch = 16, col = "black", 
      xlab = "Months", ylab = "Toxicity", 
      main = paste("Toxicity Over Time (", method, " Interpolation)", sep = "")
    )
    lines(interpolated_data$month, interpolated_data$toxicity, col = "blue", lty = 1)
    legend(
      "bottomright", 
      legend = c("Original Data", "Interpolated Data"), 
      col = c("black", "blue"), 
      pch = c(16, NA), lty = c(NA, 1)
    )
  }
  
  # Validate inputs
  validate_inputs(prevalence_data , interpolation_method)
  
  # Perform interpolation
  interpolated_data <- perform_interpolation(interpolation_method, prevalence_data , monthly_cycles)
  
  # Optionally plot results
  if (plot_output) {
    plot_results(prevalence_data , interpolated_data, interpolation_method)
  }
  
  return(interpolated_data)
}
