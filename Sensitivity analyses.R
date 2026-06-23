#### Description: Sensitivity analyses of a State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# Clear workspace
rm(list = ls())

# Load custom functions
source("Model setup.R")            # Model setup and definitions

#### Model Inputs ----
df_input_owsa <- f_input(n_sim = n_sim, setting = n_setting) # generate input parameters (min and max values for OWSA) with default n_sim
df_input_owsa <- df_input_owsa[, apply(df_input_owsa, 2, var, na.rm = TRUE) != 0] # only select parameters with variance > 0

#### Obtain deterministic sensitivity analyses ----
##### Obtain OWSA object ##### 
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

##### Inputs for TWSA ##### 
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
    twsa_8 = c("p_AI_se_arm_lymphedema", 
               "p_AI_sp_arm_lymphedema")
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
      twsa_8 = c("p_AI_se_arm_lymphedema", 
                 "p_AI_sp_arm_lymphedema")
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
        twsa_8 = c("p_AI_se_arm_lymphedema", 
                   "p_AI_sp_arm_lymphedema")
      )  
    }

##### Obtain TWSA objects ##### 
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
##### Optimal strategy plots ##### 
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

##### One way sensitivity analyses plots ##### 
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

##### Tornado plots ##### 
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

#####  TWSA plots Cost ##### 
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

##### TWSA plots QALY ##### 
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

##### TWSA plots NMB ##### 
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

