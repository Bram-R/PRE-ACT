#### Description: Value of information analyses of a State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# Clear workspace
rm(list = ls())

# Load custom functions
source("Model setup.R")            # Model setup and definitions
library("voi")                     # Load voi package

#### Run probabilistic model and obtain bcea object ----
# Create a dataframe for probabilistic sensitivity analysis (PSA) inputs
df_input <- f_input(n_sim = n_sim, setting = n_setting)

# Use a matrix to store results 
m_results_voi <- matrix(
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

# Run probabilistic analyses
for (x in 1:n_sim) m_results_voi[x, ] <- f_model(df_input[x, ])

# obtain bcea object
obj_bcea <- bcea( # create bcea object
  e = m_results_voi[,(n_treatments + 1):(n_treatments * 2)], # matrix of effects
  c = m_results_voi[,1:n_treatments], # matrix of costs
  ref = 2, # selects the 2nd row of (e,c) as reference intervention
  interventions = v_treatments, # vector of strategy names
  Kmax = 100000 # maximum value possible for the wtp
) # bcea end

##### Parameters for EVPPI ##### 
l_evppi_params <- list(
  AI = c("p_AI_se_arm_lymphedema", "p_AI_sp_arm_lymphedema", 
         "cost_t1", "cost_t2"),    
  AI_accuracy = c("p_AI_se_arm_lymphedema", "p_AI_sp_arm_lymphedema"),
  lymphedema_incidence = c("p_arm_lymphedema_m0", "p_arm_lymphedema_m2", 
                           "p_arm_lymphedema_m12", "p_arm_lymphedema_m24",
                           "p_arm_lymphedema_m36", "p_arm_lymphedema_m48", 
                           "p_arm_lymphedema_m60", "p_arm_lymphedema_m72"),
  lymphedema_consequences = c("disutility_arm_lymphedema", "cost_arm_lymphedema"),
  arm_sleeve = c("p_CP_uptake_prev_arm_lymphedema", "hr_prev_arm_lymphedema",
                 "cost_prev_arm_lymphedema_event"),
  disease_progression = c("tp_ef_ef", "p_event_lrr", "p_event_death", "tp_lrr_dm", 
                          "tp_lrr_death", "tp_dm_dm"),                                         
  costs = c("cost_t1", "cost_t2",  "cost_ef_y1_3", "cost_ef_y4_5", 
            "cost_ef_y6", "cost_lrr", "cost_dm", "cost_death", 
            "cost_lrr_event", "cost_dm_event", "cost_death_event", 
            "cost_prev_arm_lymphedema_event", "cost_prev_pain_event", 
            "cost_prev_fatigue_event", "cost_prev_fibrosis_induration_event", 
            "cost_arm_lymphedema", "cost_pain", "cost_fatigue", 
            "cost_fibrosis_induration", "cost_arm_lymphedema_event", 
            "cost_pain_event", "cost_fatigue_event", 
            "cost_fibrosis_induration_event" ),   
  utility = c("utility_ef", "utility_lrr", "utility_dm", "utility_death")            
)

#### EVPPI estimated using nonparametric regression ---- 
evppi_regression <- voi::evppi(outputs = obj_bcea, inputs = df_input, 
                               pars = l_evppi_params, check = TRUE)

# Assess regression diagnostics
check_regression(evppi_regression)

ggplot(evppi_regression, aes(x = k, y = evppi, colour = pars)) + 
  geom_line(size = 0.8) +
  labs(x = "Willingness-to-pay threshold (€)", y = "EVPPI (€)", colour = NULL) +
  theme_minimal()

### ADD ALTERNATIVE REGRESSION METHODS like "gam", "earth", "bart", "gp" to check the EVPPI estimates are similar across method

#### EVPPI estimated using traditional Monte Carlo nested loop method ---- 
f_model_voi <- f_make_model_voi()

# check with: do.call(f_model_voi, as.list(df_input[1, ]))
set.seed(12345)
evppi_mc <- voi::evppi_mc(model_fn = f_model_voi, par_fn = f_input, 
                          pars = "p_AI_se_arm_lymphedema", k = n_wtp,
                          ninner = 1000, nouter = 100)

#### EVSI estimated using nonparametric regression ---- 
evsi_regression <- voi::evsi(outputs = obj_bcea, inputs = df_input, 
                             study = "binary", 
                             n = c(50, 100, 250, 500, 1000, 2000), 
                             pars = "p_AI_se_arm_lymphedema",
                             check = TRUE)

ggplot(evsi_regression, aes(x = k, y = evsi, colour = factor(n))) +
  geom_line(size = 0.8) +
  labs(x = "Willingness-to-pay threshold (€)", y = "EVSI (€)", colour = "Sample size") +
  theme_minimal()

# Assess regression diagnostics
check_regression(evsi_regression)
