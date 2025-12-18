#### Description: Model setup, only including definitions and setup ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# General settings
options(scipen = 999, max.print = 10000, digits = 4)

# Check whether plots and text dirs are available, else create these
dirs <- c("plots", "text")
for (d in dirs) { 
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# Load and install necessary libraries
required_packages <- c(
  "docstring", "DiagrammeR", "dampack", "BCEA", "matrixStats", "summarytools"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
suppressPackageStartupMessages(lapply(required_packages, require, character.only = TRUE))

rm(required_packages, new_packages)

# Load custom functions
source("f_stm_diagram.R")          # check function using: docstring(f_stm_diagram)
source("f_input.R")                # check function using: docstring(f_input)
source("f_model.R")                # check function using: docstring(f_model)
source("f_wrapper.R")              # check function using: docstring(f_wrapper_sa) docstring(f_wrapper_intermediate)
source("f_gen_pop_utility.R")      # check function using: docstring(f_gen_pop_utility)
source("f_gen_pop_mortality.R")    # check function using: docstring(f_gen_pop_mortality)
source("f_interpolate_toxicity.R") # check function using: docstring(f_interpolate_toxicity)

#### General Setup ----
v_states <- c("Event free", "Locoregional recurrence",  # Vector of model health states
              "Distant metastasis",  "Death")          
n_states <- length(v_states)                            # Number of health states
v_treatments <- c("Current_practice",
                  "Current_practice_with_AI")           # Vector of strategy names
n_treatments <- length(v_treatments)                    # Number of treatments
v_tox <- c("arm_lymphedema", "pain", 
           "fatigue", "breast_atrophy")                 # Vector of toxicity names
n_t <- 40 * 12                                          # Model time horizon (monthly cycle)
n_sim <- 5000                                           # Number of Monte Carlo simulations
n_age_baseline <- 60                                    # Baseline age
n_p_female <- 1.00                                      # Proportion females

