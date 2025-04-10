#### Description: Sensitivity and scenario analyses of a State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# General settings
options(scipen = 999, max.print = 10000, digits = 4)

# Load and install necessary libraries
required_packages <- c(
  "dampack"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
suppressPackageStartupMessages(lapply(required_packages, require, character.only = TRUE))

# Clear workspace
rm(list = ls())

# Load custom functions
source("f_input.R")                # check function using: docstring(f_input)
source("f_model.R")                # check function using: docstring(f_model)
source("f_wrapper.R")              # check function using: docstring(f_wrapper_sa)
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
n_t <- 40 * 12                                          # Model time horizon (monthly cycle)
n_sim <- 5000                                           # Number of Monte Carlo simulations
n_age_baseline <- 60                                    # Baseline age
n_p_female <- 1.00                                      # Proportion females
n_currency <- "Pound"
n_wtp <- 30000 # willingness to pay value
v_wtp <- seq(from = 0, to = 50000, by = 2000)  # willingness to pay vector

#### Model Inputs ----
# Create a dataframe for probabilistic sensitivity analysis (PSA) inputs
m_gen_pop_utility <- f_gen_pop_utility(n_age_baseline = n_age_baseline, n_t = n_t)
m_gen_pop_mortality <- f_gen_pop_mortality(n_age_baseline = n_age_baseline, n_t = n_t, n_p_female = n_p_female)

#df_input <- f_input(n_sim = 1)
#f_wrapper_sa(params = df_input, wtp = n_wtp)

##### Obtain deterministic sensitivity analyses ----
df_input_owsa <- f_input() # generate input parameters (min and max values for OWSA) with default n_sim
df_input_owsa <- df_input_owsa[, apply(df_input_owsa, 2, var, na.rm = TRUE) != 0] # only select parameters with variance > 0

obj_owsa_dam <- run_owsa_det( # generate dampack OWSA object
  params_range = data.frame( # dataframe to be used for OWSA
    pars = names(df_input_owsa), # parameter names
    min = stack(sapply(df_input_owsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for owsa
    max = stack(sapply(df_input_owsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for owsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1), # base-case parameters (n_sim = 1 provides deterministic parameters)
  #lapply(df_input_owsa, mean), # alternative to f_input(n_sim = 1) = mean of psa inputs 
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB", "iQALY", "iCost", "iCER", "iNMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_owsa_det

df_input_twsa <- df_input_owsa[, c("disutility_arm_lymphedema", "hr_prev_arm_lymphedema")] 
obj_twsa_dam1 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

df_input_twsa <- df_input_owsa[, c("cost_arm_lymphedema", "hr_prev_arm_lymphedema")] 
obj_twsa_dam2 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

df_input_twsa <- df_input_owsa[, c("disutility_arm_lymphedema", "cost_arm_lymphedema")] 
obj_twsa_dam3 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

df_input_twsa <- df_input_owsa[, c("disutility_arm_lymphedema", "cost_t2")] 
obj_twsa_dam4 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

df_input_twsa <- df_input_owsa[, c("hr_prev_arm_lymphedema", "cost_prev_arm_lymphedema")] 
obj_twsa_dam5 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det


df_input_twsa <- df_input_owsa[, c("tp_ef_ef", "cost_arm_lymphedema")] 
obj_twsa_dam6 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

##### View deterministic sensitivity analyses results ----
owsa_opt_strat( # gives error as no parameter leads to changes in optimal strategy as they vary
  obj_owsa_dam$owsa_QALY,
  plot_const = FALSE, # TRUE = do also plot parameters that don't lead to changes in optimal strategy as they vary
  n_x_ticks = 5
) # owsa_opt_strat end

owsa_opt_strat( # gives error as no parameter leads to changes in optimal strategy as they vary
  obj_owsa_dam$owsa_Cost,
  maximize = FALSE, # need to minimize costs
  plot_const = FALSE, # TRUE = do also plot parameters that don't lead to changes in optimal strategy as they vary
  n_x_ticks = 5
) # owsa_opt_strat end

owsa_opt_strat(
  obj_owsa_dam$owsa_NMB,
  plot_const = FALSE, # TRUE = do also plot parameters that don't lead to changes in optimal strategy as they vary
  n_x_ticks = 5
) # owsa_opt_strat end

plot(
  obj_owsa_dam$owsa_QALY,
  n_y_ticks = 3, 
  n_x_ticks = 2
) # plot end

plot(
  obj_owsa_dam$owsa_Cost,
  n_y_ticks = 3, 
  n_x_ticks = 2
) # plot end

owsa_tornado(
  obj_owsa_dam$owsa_iQALY,
  #min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end

owsa_tornado(
  obj_owsa_dam$owsa_iCost,
  #min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end

owsa_tornado(
  obj_owsa_dam$owsa_iNMB,
  #min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end

owsa_tornado(
  obj_owsa_dam$owsa_iCER,
  #min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end

plot(obj_twsa_dam1$twsa_NMB)
plot(obj_twsa_dam2$twsa_NMB)
plot(obj_twsa_dam3$twsa_NMB)
plot(obj_twsa_dam4$twsa_NMB)
plot(obj_twsa_dam5$twsa_NMB)
plot(obj_twsa_dam6$twsa_NMB)

##### Scenario analyses ----
#deterministic base-case
df_input_scen <- f_input(n_sim = 1)

m_results_scen <- f_model(df_input_scen)

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) 


# 1
df_input_scen <- f_input(n_sim = 1)
df_input_scen$cost_t2 <- 10000

m_results_scen <- f_model(df_input_scen)

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

# 2
df_input_scen <- f_input(n_sim = 1)
df_input_scen$AI_se_arm_lymphedema <- 0.70

m_results_scen <- f_model(df_input_scen)

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

# 3
df_input_scen <- f_input(n_sim = 1)
df_input_scen$AI_sp_arm_lymphedema <- 0.50

m_results_scen <- f_model(df_input_scen)

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

# 4
df_input_scen <- f_input(n_sim = 1)
df_input_scen$hr_prev_arm_lymphedema <- 0.85

m_results_scen <- f_model(df_input_scen)

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

# 5
df_input_scen <- f_input(n_sim = 1)
df_input_scen$disutility_arm_lymphedema <- -0.3

m_results_scen <- f_model(df_input_scen)

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

# 6
df_input_scen <- f_input(n_sim = 1)
df_input_scen$cost_arm_lymphedema <- 5000

m_results_scen <- f_model(df_input_scen)

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end
