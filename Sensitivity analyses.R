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
df_twsa_pairs <- if(n_setting == 1) {
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
               "cost_prev_arm_lymphedema_event"),
    twsa_7 = c("p_arm_lymphedema_m72", 
               "cost_arm_lymphedema"),
    twsa_8 = c("AI_se_arm_lymphedema", 
               "AI_sp_arm_lymphedema")
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
                 "cost_prev_arm_lymphedema_event"),
      twsa_7 = c("p_arm_lymphedema_m72", 
                 "cost_arm_lymphedema"),
      twsa_8 = c("AI_se_arm_lymphedema", 
                 "AI_sp_arm_lymphedema")
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
                   "cost_prev_arm_lymphedema_event"),
        twsa_7 = c("p_arm_lymphedema_m72", 
                   "cost_arm_lymphedema"),
        twsa_8 = c("AI_se_arm_lymphedema", 
                   "AI_sp_arm_lymphedema")
      )  
    }

# Obtain TWSA objects
df_input_twsa <- df_input_owsa[, df_twsa_pairs$twsa_1] 
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

df_input_twsa <- df_input_owsa[, df_twsa_pairs$twsa_2] 
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

df_input_twsa <- df_input_owsa[, df_twsa_pairs$twsa_3] 
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

df_input_twsa <- df_input_owsa[, df_twsa_pairs$twsa_4] 
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

df_input_twsa <- df_input_owsa[, df_twsa_pairs$twsa_5] 
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

df_input_twsa <- df_input_owsa[, df_twsa_pairs$twsa_6] 
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

df_input_twsa <- df_input_owsa[, df_twsa_pairs$twsa_7] 
obj_twsa_dam_7 <- run_twsa_det( # generate dampack TWSA object
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

df_input_twsa <- df_input_owsa[, df_twsa_pairs$twsa_8] 
obj_twsa_dam_8 <- run_twsa_det( # generate dampack TWSA object
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
png(file = paste0("plots/Setting_", n_setting, "_opt_strat_", "QALY.png"), width = 500, height = 500)
owsa_opt_strat( # gives error as no parameter leads to changes in optimal strategy as they vary
  obj_owsa_dam$owsa_QALY,
  plot_const = FALSE, # TRUE = do also plot parameters that don't lead to changes in optimal strategy as they vary
  n_x_ticks = 5
) # owsa_opt_strat end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_opt_strat_", "Cost.png"), width = 500, height = 500)
owsa_opt_strat( # gives error as no parameter leads to changes in optimal strategy as they vary
  obj_owsa_dam$owsa_Cost,
  maximize = FALSE, # need to minimize costs
  plot_const = FALSE, # TRUE = do also plot parameters that don't lead to changes in optimal strategy as they vary
  n_x_ticks = 5
) # owsa_opt_strat end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_opt_strat_", "NMB.png"), width = 500, height = 500)
owsa_opt_strat(
  obj_owsa_dam$owsa_NMB,
  plot_const = FALSE, # TRUE = do also plot parameters that don't lead to changes in optimal strategy as they vary
  n_x_ticks = 5
) # owsa_opt_strat end
dev.off()

# One way sensitivity analyses plots
png(file = paste0("plots/Setting_", n_setting, "_owsa_", "QALY.png"), width = 1000, height = 700)
plot(
  obj_owsa_dam$owsa_QALY,
  n_y_ticks = 3, 
  n_x_ticks = 2
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_owsa_", "Cost.png"), width = 1000, height = 700)
plot(
  obj_owsa_dam$owsa_Cost,
  n_y_ticks = 3, 
  n_x_ticks = 2
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_owsa_", "NMB.png"), width = 1000, height = 700)
plot(
  obj_owsa_dam$owsa_NMB,
  n_y_ticks = 3, 
  n_x_ticks = 2
) # plot end
dev.off()

# Tornado plots
png(file = paste0("plots/Setting_", n_setting, "_tornado_", "iQALY.png"), width = 700, height = 500)
owsa_tornado(
  obj_owsa_dam$owsa_iQALY,
  min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_tornado_", "iCosts.png"), width = 700, height = 500)
owsa_tornado(
  obj_owsa_dam$owsa_iCost,
  min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_tornado_", "iNMB.png"), width = 700, height = 500)
owsa_tornado(
  obj_owsa_dam$owsa_iNMB,
  min_rel_diff = 0.5 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_tornado_", "iCER.png"), width = 700, height = 500)
owsa_tornado(
  obj_owsa_dam$owsa_iCER,
  min_rel_diff = 0.05 # only plot parameters that lead to a relative change in the outcome greater than or equal to this value
) # owsa_tornado end
dev.off()

# Two way sensitivity analyses strategy plots Cost
png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "1.png"), width = 700, height = 500)
plot(obj_twsa_dam_1$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "2.png"), width = 700, height = 500)
plot(obj_twsa_dam_2$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "3.png"), width = 700, height = 500)
plot(obj_twsa_dam_3$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "4.png"), width = 700, height = 500)
plot(obj_twsa_dam_4$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "5.png"), width = 700, height = 500)
plot(obj_twsa_dam_5$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "6.png"), width = 700, height = 500)
plot(obj_twsa_dam_6$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "7.png"), width = 700, height = 500)
plot(obj_twsa_dam_7$twsa_Cost)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_cost_", "8.png"), width = 700, height = 500)
plot(obj_twsa_dam_8$twsa_Cost)
dev.off()

# Two way sensitivity analyses strategy plots QALY
png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "1.png"), width = 700, height = 500)
plot(obj_twsa_dam_1$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "2.png"), width = 700, height = 500)
plot(obj_twsa_dam_2$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "3.png"), width = 700, height = 500)
plot(obj_twsa_dam_3$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "4.png"), width = 700, height = 500)
plot(obj_twsa_dam_4$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "5.png"), width = 700, height = 500)
plot(obj_twsa_dam_5$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "6.png"), width = 700, height = 500)
plot(obj_twsa_dam_6$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "7.png"), width = 700, height = 500)
plot(obj_twsa_dam_7$twsa_QALY)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_qaly_", "8.png"), width = 700, height = 500)
plot(obj_twsa_dam_8$twsa_QALY)
dev.off()

# Two way sensitivity analyses strategy plots NMB
png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "1.png"), width = 700, height = 500)
plot(obj_twsa_dam_1$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "2.png"), width = 700, height = 500)
plot(obj_twsa_dam_2$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "3.png"), width = 700, height = 500)
plot(obj_twsa_dam_3$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "4.png"), width = 700, height = 500)
plot(obj_twsa_dam_4$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "5.png"), width = 700, height = 500)
plot(obj_twsa_dam_5$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "6.png"), width = 700, height = 500)
plot(obj_twsa_dam_6$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "7.png"), width = 700, height = 500)
plot(obj_twsa_dam_7$twsa_NMB)
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_twsa_nmb_", "8.png"), width = 700, height = 500)
plot(obj_twsa_dam_8$twsa_NMB)
dev.off()

#### Obtain deterministic scenario analyses ----
#deterministic base-case
df_input_scen <- f_input(n_sim = 1, setting = n_setting)
m_results_scen_0 <- f_model(df_input_scen)

# 1 Add WGS costs + AI accuracy (based on all sources) for arm lympedema
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$cost_t2 <- df_input_scen$cost_t2 + if(n_setting == 1) {3850.71171} else 
  if(n_setting == 2) {3907.620323} else 
    if(n_setting == 3) {4744.53801}
df_input_scen$AI_se_arm_lymphedema <- if(n_setting == 1) {0.88} else 
  if(n_setting == 2) {0.88} else 
    if(n_setting == 3) {0.88}
df_input_scen$AI_sp_arm_lymphedema <- if(n_setting == 1) {0.70} else 
  if(n_setting == 2) {0.70} else 
    if(n_setting == 3) {0.70}

m_results_scen_1 <- f_model(df_input_scen)

# 2 Adjust arm lymphedema disutility to literature value
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$disutility_arm_lymphedema <- if(n_setting == 1) {-0.099} else 
  if(n_setting == 2) {-0.099} else 
    if(n_setting == 3) {-0.099}

m_results_scen_2 <- f_model(df_input_scen)

# 3 Use of Grade 0/1 arm lymphedema costs instead of Grade 3 arm lymphedema costs
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$cost_arm_lymphedema <- if(n_setting == 1) {26.163125} else 
  if(n_setting == 2) {14.51} else 
    if(n_setting == 3) {20.59284359}
df_input_scen$cost_arm_lymphedema_event <- if(n_setting == 1) {313.9575} else 
  if(n_setting == 2) {174.12} else 
    if(n_setting == 3) {247.114123}

m_results_scen_3 <- f_model(df_input_scen)

# 4 Alternative long term arm lymphedema incidence
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$p_arm_lymphedema_m72 <- if(n_setting == 1) {0.07} else 
  if(n_setting == 2) {0.07} else 
    if(n_setting == 3) {0.07}

m_results_scen_4 <- f_model(df_input_scen)

# 5 Alternative costs of arm sleeve
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$cost_prev_arm_lymphedema_event <- if(n_setting == 1) {(539.7094201 * 0.9 * 5/12)} else 
  if(n_setting == 2) {(652.4742443 * 0.9 * 5/12)} else 
    if(n_setting == 3) {(805.311062 * 0.9 * 5/12)}

m_results_scen_5 <- f_model(df_input_scen)

# 6 Alternative effectiveness of arm sleeve
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$hr_prev_arm_lymphedema <- if(n_setting == 1) {0.7} else 
  if(n_setting == 2) {0.7} else 
    if(n_setting == 3) {0.7}

m_results_scen_6 <- f_model(df_input_scen)

# 7 Assume decreased diagnostics accuracy of the PRE-ACT AI tool
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$AI_se_arm_lymphedema <- if(n_setting == 1) {0.80} else 
  if(n_setting == 2) {0.80} else 
    if(n_setting == 3) {0.80}
df_input_scen$AI_sp_arm_lymphedema <- if(n_setting == 1) {0.60} else 
  if(n_setting == 2) {0.60} else 
    if(n_setting == 3) {0.60}

m_results_scen_7 <- f_model(df_input_scen)

# 8 Assume increased diagnostics accuracy of the PRE-ACT AI tool
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$AI_se_arm_lymphedema <- if(n_setting == 1) {0.90} else 
  if(n_setting == 2) {0.90} else 
    if(n_setting == 3) {0.90}
df_input_scen$AI_sp_arm_lymphedema <- if(n_setting == 1) {0.70} else 
  if(n_setting == 2) {0.70} else 
    if(n_setting == 3) {0.70}

m_results_scen_8 <- f_model(df_input_scen)

# 9 Assumed improved quality of life and work productivity for high-risk patients because of better BMI/physical condition
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$process_utility_tp <- if(n_setting == 1) {0.01} else
  if(n_setting == 2) {0.01} else
    if(n_setting == 3) {0.01}

df_input_scen$process_utility_fp <- if(n_setting == 1) {0.01} else
  if(n_setting == 2) {0.01} else
    if(n_setting == 3) {0.01}

df_input_scen$cost_prev_arm_lymphedema_event <- if(n_setting == 1) {539.7094201 * 5/12} else
  if(n_setting == 2) {652.4742443 * 5/12} else
    if(n_setting == 3) {-765.7001469} # original cost_prev_arm_lymphedema_event adjusted for cost savings from increased productivity 

m_results_scen_9 <- f_model(df_input_scen)

# 10 Assume disutility for high-risk patients related to anxiety when being classified as high risk
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$process_utility_tp_event <- if(n_setting == 1) {-0.001} else 
  if(n_setting == 2) {-0.001} else 
    if(n_setting == 3) {-0.001}

df_input_scen$process_utility_fp_event <- if(n_setting == 1) {-0.001} else 
  if(n_setting == 2) {-0.001} else 
    if(n_setting == 3) {-0.001}

m_results_scen_10 <- f_model(df_input_scen)

# 11 Assume disutility for patients that developed arm lymphedema but that were classified as low risk, i.e. false negative classification 
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$process_utility_fn <- if(n_setting == 1) {-0.001} else 
  if(n_setting == 2) {-0.001} else 
    if(n_setting == 3) {-0.001}

m_results_scen_11 <- f_model(df_input_scen)

# 12 Add unforeseen or additional organizational and training costs 
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$cost_t2 <- df_input_scen$cost_t2 + if(n_setting == 1) {22.05370746} else 
  if(n_setting == 2) {36.2043966} else 
    if(n_setting == 3) {37.65182529}

m_results_scen_12 <- f_model(df_input_scen)

# 13 Assume utility increment for high-risk patients 
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$process_utility_tp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_scen$process_utility_fp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

m_results_scen_13 <- f_model(df_input_scen)

# 14 Assume utility increment for all patients 
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$process_utility_tp_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_scen$process_utility_fp_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_scen$process_utility_fn_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_scen$process_utility_tn_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

m_results_scen_14 <- f_model(df_input_scen)

# 15 Assume increased effectiveness of arm sleeve
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$hr_prev_arm_lymphedema <- if(n_setting == 1) {0.55} else 
  if(n_setting == 2) {0.55} else 
    if(n_setting == 3) {0.55}

m_results_scen_15 <- f_model(df_input_scen)

# 16 Combine utility beyond health scenarios 13 and 15
df_input_scen <- f_input(n_sim = 1, setting = n_setting)

df_input_scen$process_utility_tp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_scen$process_utility_fp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_scen$hr_prev_arm_lymphedema <- if(n_setting == 1) {0.55} else 
  if(n_setting == 2) {0.55} else 
    if(n_setting == 3) {0.55}

m_results_scen_16 <- f_model(df_input_scen)


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
cat("Deterministic scenario 7")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_7[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_7[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 8")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_8[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_8[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 9")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_9[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_9[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 10")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_10[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_10[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 11")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_11[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_11[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 12")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_12[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_12[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 13")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_13[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_13[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 14")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_14[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_14[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 15")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_15[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_15[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
cat("Deterministic scenario 16")
cat("\n")

calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_scen_16[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_scen_16[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

cat("\n")
sink()