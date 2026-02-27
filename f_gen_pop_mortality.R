f_gen_pop_mortality <- function(n_age_baseline, n_t, n_p_female, setting = 1) {
  #' @title Prepare Mortality Data for State Transition Model
  #' @description Processes general population mortality data to match the cycles of a 
  #' state transition model, including gender-weighted average mortality rates and their 
  #' conversion to probabilities.
  #' @param n_age_baseline Integer. The starting age of the cohort.
  #' @param n_t Integer. The total number of cycles (months) in the model.
  #' @param n_p_female Numeric. The proportion of females in the cohort (0 to 1).
  #' @param setting Integer. Selected setting `1` = UK, `2` = FR, `3` = NL. Default is `1`.
  #' @return A data frame containing filtered and processed mortality data with
  #'   columns for age and mortality probability (`prob`).
  #' @examples
  #' # Example usage:
  #' f_gen_pop_mortality(n_age_baseline = 60, n_t = 240, n_p_female = 0.5)
  
  # Load the general population mortality data
  # Data source: 2018-2020 UK life tables from ONS
  # (https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/nationallifetablesunitedkingdomreferencetables)
  mortality_data <- read.csv(
    "Gen_pop_mortality.csv",
    fileEncoding = "UTF-8-BOM"
  )

  # Calculate the gender-weighted average mortality probability per setting
  if (setting == 1) {
   
     # Calculate the gender-weighted average mortality rate
    mortality_data$mx <- mortality_data$mx_males_UK * (1 - n_p_female) + 
      mortality_data$mx_females_UK * n_p_female
    
    # Convert mortality rate to monthly probability
    mortality_data$prob <- 1 - exp(-mortality_data$mx / 12)
    
  } else if (setting == 2) {
    
    # Calculate the gender-weighted average mortality rate
    mortality_data$mx <- mortality_data$mx_males_FR * (1 - n_p_female) + 
      mortality_data$mx_females_FR * n_p_female
    
    # Convert mortality rate to monthly probability
    mortality_data$prob <- 1 - exp(-mortality_data$mx / 12)
    
  } else if (setting == 3) {
    
    # Calculate the gender-weighted average mortality probability
    mortality_data$prob <- mortality_data$p_death_males_annual_NL * (1 - n_p_female) + 
      mortality_data$p_death_females_annual_NL * n_p_female
    
    # Convert annual mortality probability to monthly probability
    mortality_data$prob <- 1 - (1 - mortality_data$prob) ^ (1/12) 
    
  } else {
    stop("Invalid setting: must be 1 (UK), 2 (FR), or 3 (NL)")
  }

  # Filter mortality data for the specified age range
  filtered_data <- mortality_data[
    mortality_data$age >= n_age_baseline & mortality_data$age <= (n_age_baseline + n_t / 12 - 1), 
    c("age", "prob")
  ]
    
  # Expand filtered data to match monthly cycles
  expanded_data <- filtered_data[
    rep(seq_len(nrow(filtered_data)), each = 12),  # Replicate each row 12 times
  ]
  
  # Add a column for the cycle number
  expanded_data$cycle <- seq_len(nrow(expanded_data))
  
  # Return the processed data
  return(expanded_data) 
}
