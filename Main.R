#### Description: Probabilistic State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# Clear workspace
rm(list = ls())

# Load custom functions
source("Model setup.R")            # Model setup and definitions

#### Model structure ----
state_labels <- list(
  state1 = c(v_states[1], "2,2!"),
  state2 = c(v_states[2], "0,0!"),
  state3 = c(v_states[3], "4,0!"),
  state4 = c(v_states[4], "2,-2!")
)

transitions <- data.frame(
  from = c("state1", "state1", "state1", "state1",
           "state2", "state2", "state2",
           "state3", "state3",
           "state4"),
  to = c("state1", "state2", "state3", "state4",
         "state2", "state3", "state4",
         "state3", "state4",
         "state4"),
  label = rep("", 10) 
)

f_stm_diagram(state_labels, transitions)
rm(state_labels, transitions)
rm(validate_state_labels, validate_transitions, generate_nodes, generate_edges, get_self_loop_style)

#### Model Inputs ----
# Country specific
n_currency <- "Pound"
n_wtp <- 30000 # willingness to pay value
v_wtp <- seq(from = 0, to = 50000, by = 2000)  # willingness to pay vector

# Create a dataframe for probabilistic sensitivity analysis (PSA) inputs
m_gen_pop_utility <- f_gen_pop_utility(n_age_baseline = n_age_baseline, n_t = n_t)
m_gen_pop_mortality <- f_gen_pop_mortality(n_age_baseline = n_age_baseline, n_t = n_t, n_p_female = n_p_female)
df_input <- f_input(n_sim = n_sim)

inputs_overview <- data.frame( # dataframe to be used for OWSA
  pars = names(df_input), # parameter names
  deterministic = as.numeric(f_input(n_sim = 1)), # deterministic values
  mean = as.numeric(lapply(df_input, mean)), # mean values
  sd = as.numeric(lapply(df_input, sd)), # sd values
  CI_min = stack(sapply(df_input, quantile, prob = 0.025, names = FALSE))[,1], # 95% percentiles 
  CI_max = stack(sapply(df_input, quantile, prob = 0.975, names = FALSE))[,1], # 95% percentiles 
  min = as.numeric(lapply(df_input, min)), # min values
  max = as.numeric(lapply(df_input, max)) # max values
) # close params_range dataframe 

sink(file = paste0("text/", "Inputs_overview.txt"))
cat("\n")
inputs_overview
cat("\n")
sink()

#### Model Results ----
# f_model(df_input[1, ]) # test model

# Use a matrix to store results 
m_results <- matrix(
  data = NA,
  nrow = n_sim,
  ncol = 3 * n_treatments,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments), paste0("LY_", v_treatments)))
)

# Run PSA
for (x in 1:n_sim) m_results[x, ] <- f_model(df_input[x, ])

# obtain dampack object
v_out_mean <- as.vector(colMeans(m_results[, 1:4])) # calculate average results

obj_icers <- calculate_icers( # create calculate_icers object
  cost = v_out_mean[1:2], # mean costs per strategy
  effect = v_out_mean[3:4], # mean effects per strategy
  strategies = v_treatments # vector of strategy names
) # calculate_icers end

obj_psa_dam <- make_psa_obj( # generate dampack PSA object
  cost = as.data.frame(m_results[,1:2]), # dataframe of costs
  effectiveness = as.data.frame(m_results[,3:4]), # dataframe of effects
  parameters = df_input, # dataframe of input parameters
  strategies = v_treatments, # vector of strategy names
) # end make_psa_obj

# obtain bcea object
obj_bcea <- bcea( # create bcea object
  e = m_results[,3:4], # matrix of effects
  c = m_results[,1:2], # matrix of costs
  ref = 2, # selects the 2nd row of (e,c) as reference intervention
  interventions = v_treatments, # vector of strategy names
  Kmax = 100000 # maximum value possible for the wtp
) # bcea end

#### Convergence statistics ----
# Running mean plots
for (x in 1:dim(m_results)[2]) {
  v_result <- m_results[, x]
  
  png(file = paste0("plots/", "convergence_", colnames(m_results)[x], ".png"), width = 1500, height = 1500)
  plot(cumsum(v_result) / seq_along(v_result), 
       type = "l",
       xlab = "Number of simulations",
       ylab = "Running mean",
       xlim = c(1000, length(v_result)),
       ylim = c(mean(v_result) - sd(v_result) * 0.5,
                mean(v_result) + sd(v_result) * 0.5),
       main = paste("Convergence plot -", colnames(m_results)[x]))
  
  abline(h = mean(v_result), col = "red", lty = 2)
  dev.off()
}

# Running incremental mean plots
m_inc_results <- matrix(
  data = cbind(
    m_results[, 1] - m_results[, 2],                   # iCosts
    m_results[, 3] - m_results[, 4],                   # iQALYs
    m_results[, 5] - m_results[, 6],                   # iLYs
    (m_results[, 1] - m_results[, 2]) /
      (m_results[, 3] - m_results[, 4])                # iCER
  ),
  nrow = n_sim,
  ncol = 4,
  dimnames = list(1:n_sim, c("iCosts", "iQALYs", "iLYs", "iCER"))
)

for (x in 1:dim(m_inc_results)[2]) {
  v_result <- m_inc_results[, x]
  
  png(file = paste0("plots/", "convergence_", colnames(m_inc_results)[x], ".png"), width = 1500, height = 1500)
  plot(cumsum(v_result) / seq_along(v_result), 
       type = "l",
       xlab = "Number of simulations",
       ylab = "Running mean",
       xlim = c(1000, length(v_result)),
       ylim = c(mean(v_result) - sd(v_result) * 0.5,
                mean(v_result) + sd(v_result) * 0.5),
       main = paste("Convergence plot -", colnames(m_inc_results)[x]))
  
  abline(h = mean(v_result), col = "red", lty = 2)
  dev.off()
}

rm(v_result, m_inc_results)

#### View results ----
summary(m_results[,1:4]) # alternatively use summarytools::dfSummary(m_results[,1:4])
plot(density(m_results[,1]))
plot(density(m_results[,2]))
plot(density(m_results[,3]))
plot(density(m_results[,4]))

sink(file = paste0("text/", "Probabilistic_base_case.txt"))
cat("\n")
obj_icers
cat("\n")
sink()

# Optionally obtain 95CI for (incremental) costs and QALYs
# costCP_CrI <- quantile(m_results[,1], probs = c(0.025, 0.975))
# costCPAI_CrI <- quantile(m_results[,2], probs = c(0.025, 0.975))
# qalyCP_CrI <- quantile(m_results[,3], probs = c(0.025, 0.975))
# qalyCPAI_CrI <- quantile(m_results[,4], probs = c(0.025, 0.975))
# icost_CrI <- quantile(m_results[,1] - m_results[,2], probs = c(0.025, 0.975))
# iqaly_CrI <- quantile(m_results[,3] - m_results[,4], probs = c(0.025, 0.975))

# Cost effectiveness frontier
png(file = paste0("plots/", "ce_frontier", ".png"), width = 1500, height = 1500)
plot(
  x = obj_icers, # icers object
  currency = n_currency, # costs units
  effect_units = "QALYs", # effects units
  label = "all" # add label to all strategies
) # plot end
dev.off()

# Incremental cost effectiveness plane
png(file = paste0("plots/", "ice_plane", ".png"), width = 1500, height = 1500)
contour2(
  he = obj_bcea, # bcea object
  wtp = n_wtp, # selects the wtp
  comparison = 1, # if more than 2 interventions, selects the pairwise comparison
  #xlim = c(-0.0005, 0.0005),
  nlevels = 4, # selects the number of levels to be plotted (default=4)
  levels = NULL, # specifies the actual levels to be plotted (default=NULL, so that R will decide)
  scale = 0.5, # scales the bandwiths for both x- and y-axis (default=0.5)
  graph = "base" # uses base graphics to produce the plot,
) # contour2 end
dev.off()

# Incremental Benefit plots
png(file = paste0("plots/", "inmb_vs_wtp", ".png"), width = 1500, height = 1500)
eib.plot(
  he = obj_bcea, # bcea object
  comparison = 1, # if more than 2 interventions, selects the pairwise comparison
  plot.cri = TRUE,
  graph = "base" # use base graphics (instead of ggplot)
) # eib.plot end
dev.off()

png(file = paste0("plots/", "inmb_density", ".png"), width = 1500, height = 1500)
ib.plot(
  he = obj_bcea, # bcea object
  comparison = 1, # if more than 2 interventions, selects the pairwise comparison
  wtp = n_wtp, # selects the wtp
  graph = "base" # uses base graphics (instead of ggplot)
) # ib.plot end
dev.off()

# CEAC
png(file = paste0("plots/", "ceac", ".png"), width = 1500, height = 1500)
plot(
  ceac(
    psa = obj_psa_dam, 
    wtp = v_wtp
  ), # ceac end
  ylim = c(0, 1),
  currency = n_currency,
  frontier = TRUE,
  points = TRUE
) # plot end
dev.off()

# ELC
png(file = paste0("plots/", "elc", ".png"), width = 1500, height = 1500)
plot(
  calc_exp_loss(
    psa = obj_psa_dam, 
    wtp = v_wtp
  ), # calc_exp_loss end
  currency = n_currency,
  n_x_ticks = 8, 
  n_y_ticks = 6
) # plot end
dev.off()

# EVPI (per patient)
png(file = paste0("plots/", "individual_evpi", ".png"), width = 1500, height = 1500)
plot(
  calc_evpi(
    psa = obj_psa_dam, 
    wtp = v_wtp,
    pop = 1
  ), # calc_evpi end
  currency = n_currency,
  n_x_ticks = 8, 
  n_y_ticks = 6
) # plot end
dev.off()

# info rank
png(file = paste0("plots/", "info_rank", ".png"), width = 1500, height = 1500)
info.rank(
  parameter = 1:ncol(df_input),
  inp = createInputs(inputs = df_input, print_is_linear_comb = FALSE),
  he = obj_bcea,
  wtp = n_wtp,
  howManyPars = 10, 
  graph = "base"
) # plot end
dev.off()

