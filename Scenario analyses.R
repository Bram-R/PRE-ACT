#### Description: Scenario analyses of a State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# Clear workspace
rm(list = ls())

# Load custom functions
source("Model setup.R")            # Model setup and definitions

#### Deterministic scenario analyses ---- 
##### 0: Base-case ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)
m_results_det_scenario_0 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_0 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_0[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_0[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) 

##### 1: Add WGS costs + AI accuracy (based on all sources) for arm lymphedema ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$cost_t2 <- df_input_det_scenario$cost_t2 + if(n_setting == 1) {3850.71171} else 
  if(n_setting == 2) {3907.620323} else 
    if(n_setting == 3) {4744.53801}
df_input_det_scenario$p_AI_se_arm_lymphedema <- if(n_setting == 1) {0.88} else 
  if(n_setting == 2) {0.88} else 
    if(n_setting == 3) {0.88}
df_input_det_scenario$p_AI_sp_arm_lymphedema <- if(n_setting == 1) {0.70} else 
  if(n_setting == 2) {0.70} else 
    if(n_setting == 3) {0.70}

m_results_det_scenario_1 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_1 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_1[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_1[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 2: Adjust arm lymphedema disutility to literature value ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$disutility_arm_lymphedema <- if(n_setting == 1) {-0.099} else 
  if(n_setting == 2) {-0.099} else 
    if(n_setting == 3) {-0.099}

m_results_det_scenario_2 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_2 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_2[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_2[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 3: Use of Grade 0/1 arm lymphedema costs instead of Grade 3 arm lymphedema costs ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$cost_arm_lymphedema <- if(n_setting == 1) {26.163125} else 
  if(n_setting == 2) {14.51} else 
    if(n_setting == 3) {20.59284359}
df_input_det_scenario$cost_arm_lymphedema_event <- if(n_setting == 1) {313.9575} else 
  if(n_setting == 2) {174.12} else 
    if(n_setting == 3) {247.114123}

m_results_det_scenario_3 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_3 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_3[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_3[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 4: Alternative long term arm lymphedema incidence ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$p_arm_lymphedema_m72 <- if(n_setting == 1) {0.07} else 
  if(n_setting == 2) {0.07} else 
    if(n_setting == 3) {0.07}

m_results_det_scenario_4 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_4 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_4[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_4[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 5: Alternative costs of arm sleeve ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$cost_prev_arm_lymphedema_event <- if(n_setting == 1) {(539.7094201 * 0.9 * 5/12)} else 
  if(n_setting == 2) {(652.4742443 * 0.9 * 5/12)} else 
    if(n_setting == 3) {(805.311062 * 0.9 * 5/12)}

m_results_det_scenario_5 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_5 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_5[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_5[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 6: Alternative effectiveness of arm sleeve ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$hr_prev_arm_lymphedema <- if(n_setting == 1) {0.7} else 
  if(n_setting == 2) {0.7} else 
    if(n_setting == 3) {0.7}

m_results_det_scenario_6 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_6 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_6[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_6[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 7: Assume decreased diagnostics accuracy of the PRE-ACT AI tool ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$p_AI_se_arm_lymphedema <- if(n_setting == 1) {0.80} else 
  if(n_setting == 2) {0.80} else 
    if(n_setting == 3) {0.80}
df_input_det_scenario$p_AI_sp_arm_lymphedema <- if(n_setting == 1) {0.60} else 
  if(n_setting == 2) {0.60} else 
    if(n_setting == 3) {0.60}

m_results_det_scenario_7 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_7 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_7[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_7[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 8: Assume increased diagnostics accuracy of the PRE-ACT AI tool ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$p_AI_se_arm_lymphedema <- if(n_setting == 1) {0.90} else 
  if(n_setting == 2) {0.90} else 
    if(n_setting == 3) {0.90}
df_input_det_scenario$p_AI_sp_arm_lymphedema <- if(n_setting == 1) {0.70} else 
  if(n_setting == 2) {0.70} else 
    if(n_setting == 3) {0.70}

m_results_det_scenario_8 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_8 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_8[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_8[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 9: Assumed improved quality of life and work productivity for high-risk patients because of better BMI/physical condition ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$process_utility_tp <- if(n_setting == 1) {0.01} else
  if(n_setting == 2) {0.01} else
    if(n_setting == 3) {0.01}

df_input_det_scenario$process_utility_fp <- if(n_setting == 1) {0.01} else
  if(n_setting == 2) {0.01} else
    if(n_setting == 3) {0.01}

df_input_det_scenario$cost_prev_arm_lymphedema_event <- if(n_setting == 1) {df_input_det_scenario$cost_prev_arm_lymphedema_event} else
  if(n_setting == 2) {df_input_det_scenario$cost_prev_arm_lymphedema_event} else
    if(n_setting == 3) {-765.7001469} # original cost_prev_arm_lymphedema_event adjusted for cost savings from increased productivity 

m_results_det_scenario_9 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_9 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_9[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_9[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 10: Assume disutility for high-risk patients related to anxiety when being classified as high risk ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$process_utility_tp_event <- if(n_setting == 1) {-0.001} else 
  if(n_setting == 2) {-0.001} else 
    if(n_setting == 3) {-0.001}

df_input_det_scenario$process_utility_fp_event <- if(n_setting == 1) {-0.001} else 
  if(n_setting == 2) {-0.001} else 
    if(n_setting == 3) {-0.001}

m_results_det_scenario_10 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_10 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_10[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_10[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 11: Assume disutility for patients that developed arm lymphedema but that were classified as low risk, i.e. false negative classification ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$process_utility_fn <- if(n_setting == 1) {-0.001} else 
  if(n_setting == 2) {-0.001} else 
    if(n_setting == 3) {-0.001}

m_results_det_scenario_11 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_11 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_11[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_11[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 12: Add unforeseen or additional organizational and training costs ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$cost_t2 <- df_input_det_scenario$cost_t2 + if(n_setting == 1) {22.05370746} else 
  if(n_setting == 2) {36.2043966} else 
    if(n_setting == 3) {37.65182529}

m_results_det_scenario_12 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_12 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_12[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_12[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 13: Assume utility increment for high-risk patients ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$process_utility_tp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_det_scenario$process_utility_fp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

m_results_det_scenario_13 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_13 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_13[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_13[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 14: Assume utility increment for all patients ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$process_utility_tp_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_det_scenario$process_utility_fp_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_det_scenario$process_utility_fn_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_det_scenario$process_utility_tn_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

m_results_det_scenario_14 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_14 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_14[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_14[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 15: Assume increased effectiveness of arm slee ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$hr_prev_arm_lymphedema <- if(n_setting == 1) {0.55} else 
  if(n_setting == 2) {0.55} else 
    if(n_setting == 3) {0.55}

m_results_det_scenario_15 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_15 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_15[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_15[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 16: Combine utility beyond health scenarios 13 and 15 ##### 
df_input_det_scenario <- f_input(n_sim = 1, setting = n_setting)

df_input_det_scenario$process_utility_tp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_det_scenario$process_utility_fp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_det_scenario$hr_prev_arm_lymphedema <- if(n_setting == 1) {0.55} else 
  if(n_setting == 2) {0.55} else 
    if(n_setting == 3) {0.55}

m_results_det_scenario_16 <- f_model(df_input_det_scenario)

obj_icers_det_scenario_16 <- calculate_icers( # create calculate_icers object
  cost = as.numeric(m_results_det_scenario_16[1:2]), # mean costs per strategy
  effect = as.numeric(m_results_det_scenario_16[3:4]), # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### Results Table (deterministic) ##### 
sink(file = paste0("text/Setting_", n_setting, "_Deterministic_scenario_analyses.txt"))
cat("\n")
cat("Deterministic base-case")
cat("\n")

obj_icers_det_scenario_0 

cat("\n")
cat("Deterministic scenario 1: Add WGS costs + AI accuracy (based on all sources) for arm lymphedema")
cat("\n")

obj_icers_det_scenario_1 

cat("\n")
cat("Deterministic scenario 2: Adjust arm lymphedema disutility to literature value")
cat("\n")

obj_icers_det_scenario_2 

cat("\n")
cat("Deterministic scenario 3: Use of Grade 0/1 arm lymphedema costs instead of Grade 3 arm lymphedema costs")
cat("\n")

obj_icers_det_scenario_3 

cat("\n")
cat("Deterministic scenario 4: Alternative long term arm lymphedema incidence")
cat("\n")

obj_icers_det_scenario_4 

cat("\n")
cat("Deterministic scenario 5: Alternative costs of arm sleeve")
cat("\n")

obj_icers_det_scenario_5 

cat("\n")
cat("Deterministic scenario 6: Alternative effectiveness of arm sleeve")
cat("\n")

obj_icers_det_scenario_6

cat("\n")
cat("Deterministic scenario 7: Assume decreased diagnostics accuracy of the PRE-ACT AI tool")
cat("\n")

obj_icers_det_scenario_7

cat("\n")

cat("\n")
cat("Deterministic scenario 8: Assume increased diagnostics accuracy of the PRE-ACT AI tool")
cat("\n")

obj_icers_det_scenario_8

cat("\n")

cat("\n")
cat("Deterministic scenario 9: Assumed improved quality of life and work productivity for high-risk patients because of better BMI/physical condition")
cat("\n")

obj_icers_det_scenario_9

cat("\n")

cat("\n")
cat("Deterministic scenario 10: Assume disutility for high-risk patients related to anxiety when being classified as high risk")
cat("\n")

obj_icers_det_scenario_10

cat("\n")

cat("\n")
cat("Deterministic scenario 11: Assume disutility for patients that developed arm lymphedema but that were classified as low risk, i.e. false negative classification")
cat("\n")

obj_icers_det_scenario_11

cat("\n")

cat("\n")
cat("Deterministic scenario 12: Add unforeseen or additional organizational and training costs")
cat("\n")

obj_icers_det_scenario_12

cat("\n")

cat("\n")
cat("Deterministic scenario 13: Assume utility increment for high-risk patients")
cat("\n")

obj_icers_det_scenario_13

cat("\n")

cat("\n")
cat("Deterministic scenario 14: Assume utility increment for all patients")
cat("\n")

obj_icers_det_scenario_14

cat("\n")

cat("\n")
cat("Deterministic scenario 15: Assume increased effectiveness of arm sleeve")
cat("\n")

obj_icers_det_scenario_15

cat("\n")

cat("\n")
cat("Deterministic scenario 16: Combine utility beyond health scenarios 13 and 15")
cat("\n")

obj_icers_det_scenario_16

cat("\n")
sink()

##### Results Figures (deterministic) #####
# Cost effectiveness frontiers for scenarios
png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_0", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_0, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_1", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_1, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_2", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_2, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_3", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_3, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_4", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_4, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_5", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_5, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_6", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_6, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_7", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_7, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_8", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_8, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_9", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_9, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_10", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_10, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_11", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_11, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_12", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_12, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_13", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_13, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_14", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_14, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_15", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_15, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_det_scenario_16", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_det_scenario_16, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

#### Probabilistic scenario analyses ---- ##CHECK WHETHER SOME PARAMETERS SHOULD BE PROBABLISTIC (e.g. using generate_gamma(4000, 500, n_sim, TRUE))
##### 0: Base-case ##### 
m_results_prob_scenario_0 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

for (x in 1:n_sim) m_results_prob_scenario_0[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_0[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_0 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 1: Add WGS costs + AI accuracy (based on all sources) for arm lymphedema ##### 
m_results_prob_scenario_1 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$cost_t2 <- df_input_prob_scenario$cost_t2 + if(n_setting == 1) {3850.71171} else 
  if(n_setting == 2) {3907.620323} else 
    if(n_setting == 3) {4744.53801}
df_input_prob_scenario$p_AI_se_arm_lymphedema <- if(n_setting == 1) {0.88} else 
  if(n_setting == 2) {0.88} else 
    if(n_setting == 3) {0.88}
df_input_prob_scenario$p_AI_sp_arm_lymphedema <- if(n_setting == 1) {0.70} else 
  if(n_setting == 2) {0.70} else 
    if(n_setting == 3) {0.70}

for (x in 1:n_sim) m_results_prob_scenario_1[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_1[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_1 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 2: Adjust arm lymphedema disutility to literature value ##### 
m_results_prob_scenario_2 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$disutility_arm_lymphedema <- if(n_setting == 1) {-0.099} else 
  if(n_setting == 2) {-0.099} else 
    if(n_setting == 3) {-0.099}

for (x in 1:n_sim) m_results_prob_scenario_2[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_2[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_2 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 3: Use of Grade 0/1 arm lymphedema costs instead of Grade 3 arm lymphedema costs ##### 
m_results_prob_scenario_3 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$cost_arm_lymphedema <- if(n_setting == 1) {26.163125} else 
  if(n_setting == 2) {14.51} else 
    if(n_setting == 3) {20.59284359}
df_input_prob_scenario$cost_arm_lymphedema_event <- if(n_setting == 1) {313.9575} else 
  if(n_setting == 2) {174.12} else 
    if(n_setting == 3) {247.114123}

for (x in 1:n_sim) m_results_prob_scenario_3[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_3[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_3 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 4: Alternative long term arm lymphedema incidence ##### 
m_results_prob_scenario_4 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$p_arm_lymphedema_m72 <- if(n_setting == 1) {0.07} else 
  if(n_setting == 2) {0.07} else 
    if(n_setting == 3) {0.07}

for (x in 1:n_sim) m_results_prob_scenario_4[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_4[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_4 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 5: Alternative costs of arm sleeve ##### 
m_results_prob_scenario_5 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$cost_prev_arm_lymphedema_event <- if(n_setting == 1) {(539.7094201 * 0.9 * 5/12)} else 
  if(n_setting == 2) {(652.4742443 * 0.9 * 5/12)} else 
    if(n_setting == 3) {(805.311062 * 0.9 * 5/12)}

for (x in 1:n_sim) m_results_prob_scenario_5[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_5[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_5 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 6: Alternative effectiveness of arm sleeve ##### 
m_results_prob_scenario_6 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$hr_prev_arm_lymphedema <- if(n_setting == 1) {0.7} else 
  if(n_setting == 2) {0.7} else 
    if(n_setting == 3) {0.7}

for (x in 1:n_sim) m_results_prob_scenario_6[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_6[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_6 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 7: Assume decreased diagnostics accuracy of the PRE-ACT AI tool ##### 
m_results_prob_scenario_7 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$p_AI_se_arm_lymphedema <- if(n_setting == 1) {0.80} else 
  if(n_setting == 2) {0.80} else 
    if(n_setting == 3) {0.80}
df_input_prob_scenario$p_AI_sp_arm_lymphedema <- if(n_setting == 1) {0.60} else 
  if(n_setting == 2) {0.60} else 
    if(n_setting == 3) {0.60}

for (x in 1:n_sim) m_results_prob_scenario_7[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_7[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_7 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 8: Assume increased diagnostics accuracy of the PRE-ACT AI tool ##### 
m_results_prob_scenario_8 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$p_AI_se_arm_lymphedema <- if(n_setting == 1) {0.90} else 
  if(n_setting == 2) {0.90} else 
    if(n_setting == 3) {0.90}
df_input_prob_scenario$p_AI_sp_arm_lymphedema <- if(n_setting == 1) {0.70} else 
  if(n_setting == 2) {0.70} else 
    if(n_setting == 3) {0.70}

for (x in 1:n_sim) m_results_prob_scenario_8[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_8[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_8 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 9: Assumed improved quality of life and work productivity for high-risk patients because of better BMI/physical condition ##### 
m_results_prob_scenario_9 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$process_utility_tp <- if(n_setting == 1) {0.01} else
  if(n_setting == 2) {0.01} else
    if(n_setting == 3) {0.01}

df_input_prob_scenario$process_utility_fp <- if(n_setting == 1) {0.01} else
  if(n_setting == 2) {0.01} else
    if(n_setting == 3) {0.01}

df_input_prob_scenario$cost_prev_arm_lymphedema_event <- if(n_setting == 1) {df_input_prob_scenario$cost_prev_arm_lymphedema_event} else
  if(n_setting == 2) {df_input_prob_scenario$cost_prev_arm_lymphedema_event} else
    if(n_setting == 3) {-765.7001469} # original cost_prev_arm_lymphedema_event adjusted for cost savings from increased productivity 

for (x in 1:n_sim) m_results_prob_scenario_9[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_9[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_9 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 10: Assume disutility for high-risk patients related to anxiety when being classified as high risk ##### 
m_results_prob_scenario_10 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$process_utility_tp_event <- if(n_setting == 1) {-0.001} else 
  if(n_setting == 2) {-0.001} else 
    if(n_setting == 3) {-0.001}

df_input_prob_scenario$process_utility_fp_event <- if(n_setting == 1) {-0.001} else 
  if(n_setting == 2) {-0.001} else 
    if(n_setting == 3) {-0.001}

for (x in 1:n_sim) m_results_prob_scenario_10[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_10[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_10 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 11: Assume disutility for patients that developed arm lymphedema but that were classified as low risk, i.e. false negative classification ##### 
m_results_prob_scenario_11 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$process_utility_fn <- if(n_setting == 1) {-0.001} else 
  if(n_setting == 2) {-0.001} else 
    if(n_setting == 3) {-0.001}

for (x in 1:n_sim) m_results_prob_scenario_11[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_11[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_11 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 12: Add unforeseen or additional organizational and training costs ##### 
m_results_prob_scenario_12 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$cost_t2 <- df_input_prob_scenario$cost_t2 + if(n_setting == 1) {22.05370746} else 
  if(n_setting == 2) {36.2043966} else 
    if(n_setting == 3) {37.65182529}

for (x in 1:n_sim) m_results_prob_scenario_12[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_12[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_12 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 13: Assume utility increment for high-risk patients ##### 
m_results_prob_scenario_13 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$process_utility_tp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_prob_scenario$process_utility_fp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

for (x in 1:n_sim) m_results_prob_scenario_13[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_13[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_13 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 14: Assume utility increment for all patients ##### 
m_results_prob_scenario_14 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$process_utility_tp_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_prob_scenario$process_utility_fp_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_prob_scenario$process_utility_fn_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_prob_scenario$process_utility_tn_event <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

for (x in 1:n_sim) m_results_prob_scenario_14[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_14[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_14 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 15: Assume increased effectiveness of arm slee ##### 
m_results_prob_scenario_15 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$hr_prev_arm_lymphedema <- if(n_setting == 1) {0.55} else 
  if(n_setting == 2) {0.55} else 
    if(n_setting == 3) {0.55}

for (x in 1:n_sim) m_results_prob_scenario_15[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_15[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_15 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### 16: Combine utility beyond health scenarios 13 and 15 ##### 
m_results_prob_scenario_16 <- matrix( # Matrix to store result 
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

df_input_prob_scenario <- f_input(n_sim = n_sim, setting = n_setting)

df_input_prob_scenario$process_utility_tp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_prob_scenario$process_utility_fp <- if(n_setting == 1) {0.001} else 
  if(n_setting == 2) {0.001} else 
    if(n_setting == 3) {0.001}

df_input_prob_scenario$hr_prev_arm_lymphedema <- if(n_setting == 1) {0.55} else 
  if(n_setting == 2) {0.55} else 
    if(n_setting == 3) {0.55}

for (x in 1:n_sim) m_results_prob_scenario_16[x, ] <- f_model(df_input_prob_scenario[x, ])

v_out_mean <- as.vector(colMeans(m_results_prob_scenario_16[, 1:(n_treatments * 2)])) # calculate average results

obj_icers_prob_scenario_16 <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:n_treatments], # mean costs per strategy
  effect = v_out_mean[(n_treatments + 1):(n_treatments * 2)], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

##### Results Table (probabilistic) ##### 
sink(file = paste0("text/Setting_", n_setting, "_Probabilistic_scenario_analyses.txt"))
cat("\n")
cat("Probabilistic base-case")
cat("\n")

obj_icers_prob_scenario_0

cat("\n")
cat("Probabilistic scenario analysis 1: Add WGS costs + AI accuracy (based on all sources) for arm lymphedema")
cat("\n")

obj_icers_prob_scenario_1

cat("\n")
cat("Probabilistic scenario analysis 2: Adjust arm lymphedema disutility to literature value")
cat("\n")

obj_icers_prob_scenario_2

cat("\n")
cat("Probabilistic scenario analysis 3: Use of Grade 0/1 arm lymphedema costs instead of Grade 3 arm lymphedema costs")
cat("\n")

obj_icers_prob_scenario_3

cat("\n")
cat("Probabilistic scenario analysis 4: Alternative long term arm lymphedema incidence")
cat("\n")

obj_icers_prob_scenario_4

cat("\n")
cat("Probabilistic scenario analysis 5: Alternative costs of arm sleeve")
cat("\n")

obj_icers_prob_scenario_5

cat("\n")
cat("Probabilistic scenario analysis 6: Alternative effectiveness of arm sleeve")
cat("\n")

obj_icers_prob_scenario_6

cat("\n")
cat("Probabilistic scenario 7: Assume decreased diagnostics accuracy of the PRE-ACT AI tool")
cat("\n")

obj_icers_prob_scenario_7

cat("\n")

cat("\n")
cat("Probabilistic scenario 8: Assume increased diagnostics accuracy of the PRE-ACT AI tool")
cat("\n")

obj_icers_prob_scenario_8

cat("\n")

cat("\n")
cat("Probabilistic scenario 9: Assumed improved quality of life and work productivity for high-risk patients because of better BMI/physical condition")
cat("\n")

obj_icers_prob_scenario_9

cat("\n")

cat("\n")
cat("Probabilistic scenario 10: Assume disutility for high-risk patients related to anxiety when being classified as high risk")
cat("\n")

obj_icers_prob_scenario_10

cat("\n")

cat("\n")
cat("Probabilistic scenario 11: Assume disutility for patients that developed arm lymphedema but that were classified as low risk, i.e. false negative classification")
cat("\n")

obj_icers_prob_scenario_11

cat("\n")

cat("\n")
cat("Probabilistic scenario 12: Add unforeseen or additional organizational and training costs")
cat("\n")

obj_icers_prob_scenario_12

cat("\n")

cat("\n")
cat("Probabilistic scenario 13: Assume utility increment for high-risk patients")
cat("\n")

obj_icers_prob_scenario_13

cat("\n")

cat("\n")
cat("Probabilistic scenario 14: Assume utility increment for all patients")
cat("\n")

obj_icers_prob_scenario_14

cat("\n")

cat("\n")
cat("Probabilistic scenario 15: Assume increased effectiveness of arm sleeve")
cat("\n")

obj_icers_prob_scenario_15

cat("\n")

cat("\n")
cat("Probabilistic scenario 16: Combine utility beyond health scenarios 13 and 15")
cat("\n")

obj_icers_prob_scenario_16

cat("\n")
sink()

##### Results Figures (probabilistic) #####
# Cost effectiveness frontiers for scenarios
png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_0", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_0, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_1", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_1, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_2", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_2, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_3", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_3, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_4", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_4, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_5", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_5, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_6", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_6, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_7", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_7, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_8", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_8, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_9", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_9, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_10", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_10, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_11", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_11, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_12", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_12, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_13", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_13, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_14", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_14, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_15", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_15, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

png(file = paste0("plots/Setting_", n_setting, "_CE_frontier_prob_scenario_16", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers_prob_scenario_16, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()
