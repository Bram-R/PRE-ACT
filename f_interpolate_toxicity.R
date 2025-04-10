f_interpolate_toxicity <- function(toxicity_data, interpolation_method = "linear", 
                                   monthly_cycles = 0:60, plot_output = FALSE) {
  
  #' Interpolate Toxicity Data and Plot Results
  #'
  #' This function interpolates toxicity data over time using either linear or spline interpolation 
  #' and optionally visualizes the results using base R's `plot()`.
  #'
  #' @param toxicity_data A data frame containing `month` and `toxicity` columns representing observed toxicity levels.
  #' @param interpolation_method A string specifying the interpolation method: "linear" or "spline". Defaults to "linear".
  #' @param monthly_cycles A numeric vector specifying the months for interpolation. Defaults to `0:60`.
  #' @param plot_output Logical; if `TRUE`, generates a plot of the interpolated results. Defaults to `TRUE`.
  #' 
  #' @return A data frame containing the interpolated toxicity values for each month in `monthly_cycles`.
  #' 
  #' @examples
  #' toxicity_data <- data.frame(
  #'   month = c(0, 12, 24, 36, 60),
  #'   toxicity = c(0.1, 0.2, 0.25, 0.3, 0.35)
  #' )
  #' 
  #' result <- interpolate_toxicity(toxicity_data, interpolation_method = "spline", plot_output = TRUE)
  #' print(result)
  #'
  
  # Helper function: Validate inputs
  validate_inputs <- function(toxicity_data, interpolation_method) {
    if (!is.data.frame(toxicity_data) || !all(c("month", "toxicity") %in% names(toxicity_data))) {
      stop("`toxicity_data` must be a data frame with columns `month` and `toxicity`.")
    }
    if (!interpolation_method %in% c("linear", "spline")) {
      stop("`interpolation_method` must be either 'linear' or 'spline'.")
    }
  }
  
  # Helper function: Perform interpolation
  perform_interpolation <- function(method, toxicity_data, monthly_cycles) {
    if (method == "linear") {
      interpolated <- approx(
        x = toxicity_data$month, 
        y = toxicity_data$toxicity, 
        xout = monthly_cycles
      )
    } else if (method == "spline") {
      interpolated <- spline(
        x = toxicity_data$month, 
        y = toxicity_data$toxicity, 
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
  validate_inputs(toxicity_data, interpolation_method)
  
  # Perform interpolation
  interpolated_data <- perform_interpolation(interpolation_method, toxicity_data, monthly_cycles)
  
  # Optionally plot results
  if (plot_output) {
    plot_results(toxicity_data, interpolated_data, interpolation_method)
  }
  
  return(interpolated_data)
}
