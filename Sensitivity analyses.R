#### Description: Sensitivity and scenario analyses of a State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# Clear workspace
rm(list = ls())

# Load custom functions
source("Model setup.R")            # Model setup and definitions

#### Model Inputs ----
df_input_owsa <- f_input(n_sim = n_sim, setting = n_setting) # generate input parameters (min and max values for OWSA) with default n_sim
df_input_owsa <- df_input_owsa[, apply(df_input_owsa, 2, var, na.rm = TRUE) != 0] # only select parameters with variance > 0

#### Obtain deterministic sensitivity analyses ----
# Obtain OWSA object
obj_owsa_dam <- run_owsa_det( # generate dampack OWSA object
  params_range = data.frame( # dataframe to be used for OWSA
    pars = names(df_input_owsa), # parameter names
    min = stack(sapply(df_input_owsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for owsa
    max = stack(sapply(df_input_owsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for owsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1, setting = n_setting), # base-case parameters (n_sim = 1 provides deterministic parameters)
  #lapply(df_input_owsa, mean), # alternative to f_input(n_sim = 1, setting = n_setting) = mean of psa inputs 
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB", "iQALY", "iCost", "iCER", "iNMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_owsa_det

# Inputs for TWSA
df_input_twsa <- if(n_setting == 1) {
  data.frame(
    twsa_1 = c("cost_t2", 
               "disutility_arm_lymphedema"),
    twsa_2 = c("cost_t2", 
               "hr_prev_arm_lymphedema"),
    twsa_3 = c("cost_t2", 
               "cost_prev_arm_lymphedema_event"),
    twsa_4 = c("disutility_arm_lymphedema", 
               "hr_prev_arm_lymphedema"),
    twsa_5 = c("disutility_arm_lymphedema", 
               "cost_prev_arm_lymphedema_event"),
    twsa_6 = c("hr_prev_arm_lymphedema", 
               "cost_prev_arm_lymphedema_event")
  )} else if(n_setting == 2) {
    data.frame(
      twsa_1 = c("cost_t2", 
                 "disutility_arm_lymphedema"),
      twsa_2 = c("cost_t2", 
                 "hr_prev_arm_lymphedema"),
      twsa_3 = c("cost_t2", 
                 "cost_prev_arm_lymphedema_event"),
      twsa_4 = c("disutility_arm_lymphedema", 
                 "hr_prev_arm_lymphedema"),
      twsa_5 = c("disutility_arm_lymphedema", 
                 "cost_prev_arm_lymphedema_event"),
      twsa_6 = c("hr_prev_arm_lymphedema", 
                 "cost_prev_arm_lymphedema_event")
    )} else if(n_setting == 3) { 
      data.frame(
        twsa_1 = c("cost_t2", 
                   "disutility_arm_lymphedema"),
        twsa_2 = c("cost_t2", 
                   "hr_prev_arm_lymphedema"),
        twsa_3 = c("cost_t2", 
                   "cost_prev_arm_lymphedema_event"),
        twsa_4 = c("disutility_arm_lymphedema", 
                   "hr_prev_arm_lymphedema"),
        twsa_5 = c("disutility_arm_lymphedema", 
                   "cost_prev_arm_lymphedema_event"),
        twsa_6 = c("hr_prev_arm_lymphedema", 
                   "cost_prev_arm_lymphedema_event")
      )  
    }

# Obtain TWSA objects
df_input_twsa <- df_input_owsa[, df_input_twsa$twsa_1] 
obj_twsa_dam_1 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1, setting = n_setting), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

df_input_twsa <- df_input_owsa[, df_input_twsa$twsa_2] 
obj_twsa_dam_2 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1, setting = n_setting), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

df_input_twsa <- df_input_owsa[, df_input_twsa$twsa_3] 
obj_twsa_dam_3 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1, setting = n_setting), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

df_input_twsa <- df_input_owsa[, df_input_twsa$twsa_4] 
obj_twsa_dam_4 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1, setting = n_setting), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

df_input_twsa <- df_input_owsa[, df_input_twsa$twsa_5] 
obj_twsa_dam_5 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1, setting = n_setting), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

df_input_twsa <- df_input_owsa[, df_input_twsa$twsa_6] 
obj_twsa_dam_6 <- run_twsa_det( # generate dampack TWSA object
  params_range = data.frame( # dataframe to be used for TWSA
    pars = names(df_input_twsa), # parameter names
    min = stack(sapply(df_input_twsa, quantile, prob = 0.025, names = FALSE))[,1], # use 95% percentiles for twsa
    max = stack(sapply(df_input_twsa, quantile, prob = 0.975, names = FALSE))[,1] # use 95% percentiles for twsa 
  ), # close params_range dataframe 
  params_basecase = f_input(n_sim = 1, setting = n_setting), # base-case parameters (n_sim = 1 provides deterministic parameters)
  nsamp = 100, # number of sets of parameter values to be generated (between min and max)
  FUN = f_wrapper_sa, wtp = n_wtp, # function that produces outcomes
  outcomes = c("QALY", "Cost", "NMB"), # outcomes of interest produced by FUN 
  progress = TRUE # show progression in console
) # end run_twsa_det

#### Deterministic sensitivity analyses results ----
# Optimal strategy plots
png(file = paste0("plots/Setting_", n_setting, "_opt_strat_", "QALY.png"), width = 1500, height = 1500)
owsa_opt_strat( # gives error as no parameter leads to changes in optimal strategy as they vary
  obj_owsa_dam$owsa_QALY,
  plot_const = FALSE, # TRUE = do also plot parameters that don't lead to changes in optimal strategy as they vary
  n_x_ticks = 5
) # owsa_opt_strat end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_opt_strat_", "Cost.png"), width = 1500, height = 1500)
owsa_opt_strat( # gives error as no parameter leads to changes in optimal strategy as they vary
  obj_owsa_dam$owsa_Cost,
  maximize = FALSE, # need to minimize costs
  plot_const = FALSE, # TRUE = do also plot parameters that don't lead to changes in optimal strategy as they vary
  n_x_ticks = 5
) # owsa_opt_strat end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_opt_strat_", "NMB.png"), width = 1500, height = 1500)
owsa_opt_strat(
  obj_owsa_dam$owsa_NMB,
  plot_const = FALSE, # TRUE = do also plot parameters that don't lead to changes in optimal strategy as they vary
  n_x_ticks = 5
) # owsa_opt_strat end
dev.off()

# One way sensitivity analyses plots
png(file = paste0("plots/Setting_", n_setting, "_owsa_", "QALY.png"), width = 1500, height = 1500)
plot(
  obj_owsa_dam$owsa_QALY,
  n_y_ticks = 3, 
  n_x_ticks = 2
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_owsa_", "Cost.png"), width = 1500, height = 1500)
plot(
  obj_owsa_dam$owsa_Cost,
  n_y_ticks = 3, 
  n_x_ticks = 2
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_owsa_", "NMB.png"), width = 1500, height = 1500)
plot(
  obj_owsa_dam$owsa_NMB,
  n_y_ticks = 3, 
  n_x_ticks = 2
) # plot end
dev.off()

# Tornado plots
png(file = paste0("plots/Setting_", n_setting, "_tornado_", "iQALY.png"), width = 1500, height = 1500)
owsa_tornado(
  obj_owsa_dam$owsa_iQALY,
  #min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_tornado_", "iCosts.png"), width = 1500, height = 1500)
owsa_tornado(
  obj_owsa_dam$owsa_iCost,
  #min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_tornado_", "iNMB.png"), width = 1500, height = 1500)
owsa_tornado(
  obj_owsa_dam$owsa_iNMB,
  #min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_tornado_", "iCER.png"), width = 1500, height = 1500)
owsa_tornado(
  obj_owsa_dam$owsa_iCER,
  #min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end
dev.off()

# Two way sensitivity analyses strategy plots Cost
png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "1.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_1$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "2.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_2$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "3.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_3$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "4.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_4$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "5.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_5$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "6.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_6$twsa_Cost)
dev.off()

# Two way sensitivity analyses strategy plots QALY
png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "1.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_1$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "2.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_2$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "3.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_3$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "4.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_4$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "5.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_5$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "6.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_6$twsa_QALY)
dev.off()

# Two way sensitivity analyses strategy plots NMB
png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "1.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_1$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "2.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_2$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "3.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_3$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "4.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_4$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "5.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_5$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "6.png"), width = 1500, height = 1500)
plot(obj_twsa_dam_6$twsa_NMB)
dev.off()

#### Obtain deterministic scenario analyses ----
#deterministic base-case
df_input_scen <- f_input(n_sim = 1, setting = n_setting)
m_results_scen_0 <- f_model(df_input_scen)

# 1 Add WGS costs + AI accuracy (based on all sources) for arm lympedema
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$cost_t2 <- if(n_setting == 1) {1} else 
  if(n_setting == 2) {2} else 
    if(n_setting == 3) {3}
df_input_scen$AI_se_arm_lymphedema <- if(n_setting == 1) {1} else 
  if(n_setting == 2) {2} else 
    if(n_setting == 3) {3}
df_input_scen$AI_sp_arm_lymphedema <- if(n_setting == 1) {1} else 
  if(n_setting == 2) {2} else 
    if(n_setting == 3) {3}

m_results_scen_1 <- f_model(df_input_scen)

# 2 Adjust arm lymphedema disutility to literature value
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$disutility_arm_lymphedema <- if(n_setting == 1) {1} else 
  if(n_setting == 2) {2} else 
    if(n_setting == 3) {3}

m_results_scen_2 <- f_model(df_input_scen)

# 3 Exclusion of Grade 3 arm lymphedema costs
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$cost_arm_lymphedema <- if(n_setting == 1) {1} else 
  if(n_setting == 2) {2} else 
    if(n_setting == 3) {3}
df_input_scen$cost_arm_lymphedema_event <- if(n_setting == 1) {1} else 
  if(n_setting == 2) {2} else 
    if(n_setting == 3) {3}

m_results_scen_3 <- f_model(df_input_scen)

# 4 Alternative long term arm lymphedema incidence
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$p_arm_lymphedema_m72 <- if(n_setting == 1) {1} else 
  if(n_setting == 2) {2} else 
    if(n_setting == 3) {3}

m_results_scen_4 <- f_model(df_input_scen)

# 5 Alternative costs of arm sleeve
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$cost_prev_arm_lymphedema_event <- if(n_setting == 1) {1} else 
  if(n_setting == 2) {2} else 
    if(n_setting == 3) {3}

m_results_scen_5 <- f_model(df_input_scen)

# 6 Alternative effectiveness of arm sleeve
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$hr_prev_arm_lymphedema <- if(n_setting == 1) {1} else 
  if(n_setting == 2) {2} else 
    if(n_setting == 3) {3}

m_results_scen_6 <- f_model(df_input_scen)

#### Deterministic scenario analyses results ----
sink(file = paste0("text/Setting_", n_setting, "_Deterministic_scenario_analyses.txt"))
cat("\n")
cat("Deterministic base-case")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_0[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_0[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) 

cat("\n")
cat("Deterministic scenario 1")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_1[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_1[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 2")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_2[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_2[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 3")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_3[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_3[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 4")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_4[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_4[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 5")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_5[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_5[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 6")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_6[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_6[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
sink()
