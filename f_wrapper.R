f_wrapper_sa <- function(params, wtp = 30000) {
  #' Wrapper Function for Sensitivity Analysis
  #'
  #' This function serves as a wrapper for conducting sensitivity analyses using the `dampack` package. 
  #' It runs the `f_model()` function with a given set of parameters and structures the results 
  #' for cost-effectiveness analysis, including incremental outcomes.
  #'
  #' @param params A named list containing model parameters.
  #' @param wtp A numeric value representing the willingness-to-pay (WTP) threshold 
  #' for cost-effectiveness analysis (default: 30000).
  #'
  #' @details
  #' This function executes a single simulation of `f_model()` and extracts key cost-effectiveness 
  #' outcomes, including costs, QALYs, and incremental cost-effectiveness ratios (iCERs). 
  #' It is designed for use with deterministic one-way sensitivity analysis (OWSA) using 
  #' `dampack::run_owsa_det()`.
  #'
  #' The output includes:
  #' \itemize{
  #'   \item **Total costs and QALYs** for each treatment strategy.
  #'   \item **Net monetary benefit (NMB)** for each strategy, computed as:  
  #'         \eqn{NMB = (QALYs \times WTP) - Costs}.
  #'   \item **Incremental outcomes**, including:
  #'     \itemize{
  #'       \item Incremental QALYs (\code{iQALY}).
  #'       \item Incremental costs (\code{iCost}).
  #'       \item Incremental cost-effectiveness ratio (\code{iCER}).
  #'       \item Incremental net monetary benefit (\code{iNMB}).
  #'     }
  #' }
  #'
  #' @return A data frame with cost-effectiveness results, including:
  #' \itemize{
  #'   \item \code{Option}: Treatment strategy.
  #'   \item \code{QALY}: Quality-adjusted life years.
  #'   \item \code{Cost}: Total cost.
  #'   \item \code{NMB}: Net monetary benefit.
  #'   \item \code{iQALY}, \code{iCost}, \code{iCER}, \code{iNMB}: Incremental outcomes.
  #' }
  #'
  #' @export
  
  # Run model simulation
  v_results <- f_model(params)
  
  # Create results dataframe for cost-effectiveness analysis
  df_results <- data.frame(
    Option = v_treatments,
    QALY = c(v_results[3], v_results[4]),
    Cost = c(v_results[1], v_results[2]),
    NMB = c(v_results[3] * wtp - v_results[1], v_results[4] * wtp - v_results[2]),
    iQALY = rep(v_results[4] - v_results[3], n_treatments),
    iCost = rep(v_results[2] - v_results[1], n_treatments),
    iCER = rep((v_results[2] - v_results[1]) / (v_results[4] - v_results[3]), n_treatments),
    iNMB = rep((v_results[4] * wtp - v_results[2]) - (v_results[3] * wtp - v_results[1]), n_treatments)
  )
  
  return(df_results)
}

f_wrapper_intermediate <- function(df_input) {
  #' Wrapper Function for Extracting Intermediate Simulation Results
  #'
  #' This function runs the `f_model()` function for multiple parameter sets (simulations) 
  #' and extracts the intermediate cycle-level results into a structured 3D array.
  #'
  #' @param df_input A data frame where each row represents a set of model parameters 
  #' to be used in separate simulations.
  #'
  #' @details 
  #' This function iterates over each row in `df_input`, passing the corresponding parameters 
  #' to `f_model()` with `intermediate = TRUE`. The results from each simulation are stored 
  #' in a 3D array with the following structure:
  #' \itemize{
  #'   \item **Dimension 1 (rows)**: Simulation cycles (e.g., months or years).
  #'   \item **Dimension 2 (columns)**: Outcomes tracked, including costs, QALYs, and Markov traces/LYs.
  #'   \item **Dimension 3 (depth)**: Individual simulations (one per row in `df_input`).
  #' }
  #'
  #' The column names in the output array include:
  #' \itemize{
  #'   \item Costs and QALYs per health state for each treatment strategy.
  #'   \item Toxicity-related costs and QALYs for each treatment.
  #'   \item Event-related costs and QALYs.
  #'   \item Markov trace/LYs (state occupancy probabilities).
  #'   \item Toxicity incidence rates for each treatment.
  #' }
  #'
  #' @return A 3D array of dimension `(cycles, outcomes, simulations)`, where:
  #' \itemize{
  #'   \item **Rows** represent model cycles.
  #'   \item **Columns** represent different outcomes (costs, QALYs, Markov trace/LYs).
  #'   \item **Depth** represents different simulations (one per row in `df_input`).
  #' }
  #'
  #' @export
  
  # create array to store results; 3 dimensions (cycles, outcomes, simulations)
  a_results <- array(data = NA, 
                     dim = c(n_t + 1, # number of cycles
                             3 * n_treatments * (length(v_states) + length(v_tox) + 1), # number of results (number of outcomes * number of treatments * number of health states/events/toxicities)
                             dim(df_input)[1] # number of simulations
                     ), # close dim
                     dimnames = list( # name dimensions
                       paste0("cycle_", 1:(n_t + 1)),
                       c(paste0("Cost_", v_states, "_", v_treatments[1]),   # Cost per health state for t1
                         paste0("Cost_", v_tox, "_", v_treatments[1]),      # Tox cost for t1
                         paste0("Cost_Event_", v_treatments[1]),            # Event cost for t1
                         
                         paste0("Cost_", v_states, "_", v_treatments[2]),   # Cost per health state for t2
                         paste0("Cost_", v_tox, "_", v_treatments[2]),      # Tox cost for t2
                         paste0("Cost_Event_", v_treatments[2]),            # Event cost for t2
                         
                         paste0("QALY_", v_states, "_", v_treatments[1]),   # QALYs per health state for t1
                         paste0("QALY_", v_tox, "_", v_treatments[1]),      # Tox QALY for t1
                         paste0("QALY_Event_", v_treatments[1]),            # Event QALY for t1
                         
                         paste0("QALY_", v_states, "_", v_treatments[2]),   # QALYs per health state for t2
                         paste0("QALY_", v_tox, "_", v_treatments[2]),      # Tox QALY for t2
                         paste0("QALY_Event_", v_treatments[2]),            # Event QALY for t2
                         
                         paste0("LY_", v_states, "_", v_treatments[1]),     # LYs per health state for t1 (Markov trace)
                         paste0("Incidence_", v_tox, "_", v_treatments[1]), # Tox incidence for t1
                         paste0("LY_Event_", v_treatments[1]),              # Event LYs for t1
                         
                         paste0("LY_", v_states, "_", v_treatments[2]),     # LYs per health state for t2 (Markov trace)
                         paste0("Incidence_", v_tox, "_", v_treatments[2]), # Tox incidence for t2
                         paste0("LY_Event_", v_treatments[2])               # Event LYs for t2
                       ), # end c
                       paste0("sim_", 1:dim(df_input)[1])
                     ) # close dimnames 
  ) # end array 
  
  # Run model simulation
  for (x in 1:dim(df_input)[1]) a_results[,,x ] <- f_model(df_input[x, ], intermediate = TRUE)
  
  return(a_results)
}


