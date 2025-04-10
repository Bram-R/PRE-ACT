#### Description: Intermediate results of a State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# General settings
options(scipen = 999, max.print = 10000, digits = 4)

# Load and install necessary libraries
required_packages <- c(
  "dampack", "matrixStats", "summarytools"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
suppressPackageStartupMessages(lapply(required_packages, require, character.only = TRUE))

# Clear workspace
rm(list = ls())

# Load custom functions
source("f_input.R")                # check function using: docstring(f_input)
source("f_model.R")                # check function using: docstring(f_model)
source("f_wrapper.R")              # check function using: docstring(f_wrapper_intermediate)
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
df_input <- f_input(n_sim = n_sim)

##### Obtain intermediate results using wrapper intermediate function ----
# Probabilistic results (intermediate)
a_out_interm <- f_wrapper_intermediate(df_input)

# str(a_out_interm[,,]) # inspect
# names(a_out_interm[1, , 1]) # inspect

##### Costs and QALYs per cycle ----
a_out_cycle_res <- array(
  c(apply(a_out_interm[, 1:9, ], c(1, 3), sum),
    apply(a_out_interm[, 10:18, ], c(1, 3), sum),
    apply(a_out_interm[, 19:27, ], c(1, 3), sum),
    apply(a_out_interm[, 28:36, ], c(1, 3), sum)
  ),
  dim = c(481, 5000, 4)
)

a_out_cycle_res <- aperm(a_out_cycle_res, perm = c(1, 3, 2))

m_out_cycle_res <- rowMeans(a_out_cycle_res[, 1:4, ], dims = 2)
m_out_cycle_res_percentiles <- cbind(
  rowQuantiles(a_out_cycle_res[, 1, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_cycle_res[, 2, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_cycle_res[, 3, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_cycle_res[, 4, ], probs = c(0.025, 0.975))
)

# Probabilistic costs
matplot(
  x = 0:n_t,
  y = cbind(
    cumsum(m_out_cycle_res[, 1]),
    cumsum(m_out_cycle_res_percentiles[, 1]), 
    cumsum(m_out_cycle_res_percentiles[, 2])
  ),
  type = "l",
  ylab = "Costs",
  xlab = "Cycle",
  main = paste0("Cumulative costs (probabilistic) ", v_treatments[1]),
  col = c(1, 1, 1),
  lty = c(1, 2, 2)
)
legend(
  "topleft", 
  inset = c(0.05, 0),
  c("Probabilistic mean", dimnames(m_out_cycle_res_percentiles[, 1:2])[[2]]),
  cex = 0.5,
  col = c(1, 1, 1),
  lty = c(1, 2, 2),
  bty = "n"
)

matplot(
  x = 0:n_t,
  y = cbind(
    cumsum(m_out_cycle_res[, 2]),
    cumsum(m_out_cycle_res_percentiles[, 3]), 
    cumsum(m_out_cycle_res_percentiles[, 4])
  ),
  type = "l",
  ylab = "Costs",
  xlab = "Cycle",
  main = paste0("Cumulative costs (probabilistic) ", v_treatments[2]),
  col = c(1, 1, 1),
  lty = c(1, 2, 2)
)
legend(
  "topleft", 
  inset = c(0.05, 0),
  c("Probabilistic mean", dimnames(m_out_cycle_res_percentiles[, 1:2])[[2]]),
  cex = 0.5,
  col = c(1, 1, 1),
  lty = c(1, 2, 2),
  bty = "n"
)

matplot(
  x = 0:n_t,
  y = cbind(
    cumsum(m_out_cycle_res[, 1]),
    cumsum(m_out_cycle_res[, 2])
  ),
  type = "l",
  ylab = "Costs",
  xlab = "Cycle",
  main = paste0("Cumulative costs (probabilistic) ", v_treatments[1], " and ", v_treatments[2]),
  col = c(2, 2),
  lty = c(1, 2)
)
legend(
  "topleft", 
  inset = c(0.05, 0),
  v_treatments,
  cex = 0.5,
  col = c(2, 2),
  lty = c(1, 2),
  bty = "n"
)

# Probabilistic QALYs 
matplot(
  x = 0:n_t,
  y = cbind(
    cumsum(m_out_cycle_res[, 3]),
    cumsum(m_out_cycle_res_percentiles[, 5]), 
    cumsum(m_out_cycle_res_percentiles[, 6])
  ),
  type = "l",
  ylab = "QALYs",
  xlab = "Cycle",
  main = paste0("Cumulative QALYs (probabilistic) ", v_treatments[1]),
  col = c(1, 1, 1),
  lty = c(1, 2, 2)
)
legend(
  "bottomright", 
  inset = c(0.05, 0),
  c("Probabilistic mean", dimnames(m_out_cycle_res_percentiles[, 1:2])[[2]]),
  cex = 0.5,
  col = c(1, 1, 1),
  lty = c(1, 2, 2),
  bty = "n"
)

matplot(
  x = 0:n_t,
  y = cbind(
    cumsum(m_out_cycle_res[, 4]),
    cumsum(m_out_cycle_res_percentiles[, 7]), 
    cumsum(m_out_cycle_res_percentiles[, 8])
  ),
  type = "l",
  ylab = "QALYs",
  xlab = "Cycle",
  main = paste0("Cumulative QALYs (probabilistic) ", v_treatments[2]),
  col = c(1, 1, 1),
  lty = c(1, 2, 2)
)
legend(
  "bottomright", 
  inset = c(0.05, 0),
  c("Probabilistic mean", dimnames(m_out_cycle_res_percentiles[, 1:2])[[2]]),
  cex = 0.5,
  col = c(1, 1, 1),
  lty = c(1, 2, 2),
  bty = "n"
)

matplot(
  x = 0:n_t,
  y = cbind(
    cumsum(m_out_cycle_res[, 3]),
    cumsum(m_out_cycle_res[, 4])
  ),
  type = "l",
  ylab = "QALYs",
  xlab = "Cycle",
  main = paste0("Cumulative QALYs (probabilistic) ", v_treatments[1], " and ", v_treatments[2]),
  col = c(2, 2),
  lty = c(1, 2)
)
legend(
  "topleft", 
  inset = c(0.05, 0),
  v_treatments,
  cex = 0.5,
  col = c(2, 2),
  lty = c(1, 2),
  bty = "n"
)

##### Toxicity per cycle ----
m_toxicity <- rowMeans(a_out_interm[, 45:52, ], dims = 2)
m_toxicity_percentiles <- cbind(
  rowQuantiles(a_out_interm[, 45, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 46, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 47, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 48, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 49, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 50, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 51, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 52, ], probs = c(0.025, 0.975))
)

# Toxicity incidence
matplot(
  x = 0:n_t,
  y = m_toxicity[, 1:8], 
  ylim = c(0, 0.4),
  type = "l",
  ylab = "Toxicity incidence (conditional on being alive)",
  xlab = "Cycle",
  main = paste0("Average toxicity incidence for ", v_treatments[1], " and ", v_treatments[2]),
  col = rep(2:5, 2),
  lty = c(rep(1, 4), rep(2, 4))
) # matplot end
legend(
  "topleft", 
  inset = c(0.05, 0),
  dimnames(m_toxicity[, 1:8])[[2]], 
  cex = 0.75,
  col = rep(2:5,2), 
  lty =  c(rep(1, 4), rep(2, 4)),
  bty = "n"
) # legend end

# Arm lymphedema incidence
matplot(
  x = 0:n_t,
  y = cbind(
    m_toxicity[, 1],
    m_toxicity_percentiles[, 1], 
    m_toxicity_percentiles[, 2],
    m_toxicity[, 5],
    m_toxicity_percentiles[, 9], 
    m_toxicity_percentiles[, 10]
  ),
  type = "l",
  ylab = "Arm lymphedema incidence (conditional on being alive)",
  xlab = "Cycle",
  xlim = c(0, 100),
  main = "Arm lymphedema incidence (conditional on being alive)",
  col = c(2, 2, 2, 3, 3, 3),
  lty = c(1, 2, 2, 1, 2, 2)
)
legend(
  "topright", 
  inset = c(0.05, 0),
  c(paste0("Probabilistic mean ", v_treatments[1]), dimnames(m_out_cycle_res_percentiles[, 1:2])[[2]],
    paste0("Probabilistic mean ", v_treatments[2]), dimnames(m_out_cycle_res_percentiles[, 1:2])[[2]]),
  cex = 0.5,
  col = c(2, 2, 2, 3, 3, 3),
  lty = c(1, 2, 2, 1, 2, 2),
  bty = "n"
)

##### Traces per cycle ----
m_trace <- rowMeans(a_out_interm[, 37:44, ], dims = 2)

# Traces
matplot(
  x = 0:n_t,
  y = m_trace[, 1:8], 
  ylim = c(0, 1),
  type = "l",
  ylab = "Probability of state occupancy",
  xlab = "Cycle",
  main = paste0("Average traces for ", v_treatments[1], " and ", v_treatments[2]),
  col = rep(2:5, 2),
  lty = c(rep(1, 4), rep(2, 4))
) # matplot end
legend(
  "topleft", 
  inset = c(0.05, 0),
  dimnames(m_trace[, 1:8])[[2]], 
  cex = 0.75,
  col = rep(2:5,2), 
  lty =  c(rep(1, 4), rep(2, 4)),
  bty = "n"
) # legend end

# Proportion death
matplot(
  x = 0:n_t,
  y = cbind(m_trace[, 4], m_trace[, 8]), 
  ylim = c(0, 1),
  type = "l",
  ylab = "Proportion death",
  xlab = "Cycle",
  main = paste0("Average proportion death for ", v_treatments[1], " and ", v_treatments[2]),
  col = rep(1, 2),
  lty = c(1, 2)
) # matplot end
legend(
  "topleft", 
  inset = c(0.05, 0),
  c(names(m_trace[1,])[4], names(m_trace[1,])[8]), 
  cex = 0.75,
  col = rep(1, 2),
  lty = c(1, 2),
  bty = "n"
) # legend end

# Traces conditional on being alive
matplot(
  x = 0:n_t,
  y = cbind(m_trace[, 1:3]/(1 - m_trace[, 4]), m_trace[, 5:7]/(1 - m_trace[, 8])), 
  ylim = c(0, 1),
  type = "l",
  ylab = "Proportion (conditional on being alive)",
  xlab = "Cycle",
  main = paste0("Average traces conditional on being alive ", v_treatments[1], " and ", v_treatments[2]),
  col = rep(2:4, 2),
  lty = c(rep(1, 3), rep(2, 3))
) # matplot end
legend(
  "topleft", 
  inset = c(0.05, 0),
  c(names(m_trace[1,])[1:3], names(m_trace[1,])[5:7]), 
  cex = 0.75,
  col = rep(2:4, 2),
  lty = c(rep(1, 3), rep(2, 3)),
  bty = "n"
) # legend end

##### Disaggregated costs and QALYs per cycle ----
m_out_dis_cycle_res <- rowMeans(a_out_interm[, 1:36, ], dims = 2)

# Health state costs
matplot(
  x = 0:n_t,
  y = m_out_dis_cycle_res[, 1:18], 
  type = "l",
  ylab = "Costs",
  xlab = "Cycle",
  main = paste0("Average health state costs for ", v_treatments[1], " and ", v_treatments[2]),
  col = rep(2:4, 2),
  lty = c(rep(1, 3), rep(2, 3))
) # matplot end
legend(
  "topright", 
  inset = c(0.05, 0),
  dimnames(m_out_dis_cycle_res[, 1:18])[[2]], 
  cex = 0.75,
  col = rep(2:4,2), 
  lty =  c(rep(1, 3), rep(2, 3)),
  bty = "n"
) # legend end

# Health state QALYs
matplot(
  x = 0:n_t,
  y = m_out_dis_cycle_res[, 19:36], 
  ylim = c(0, 1),
  type = "l",
  ylab = "Costs",
  xlab = "Cycle",
  main = paste0("Average health state QALYs ", v_treatments[1], " and ", v_treatments[2]),
  col = rep(2:4, 2),
  lty = c(rep(1, 3), rep(2, 3))
) # matplot end
legend(
  "topright", 
  inset = c(0.05, 0),
  dimnames(m_out_dis_cycle_res[, 19:36])[[2]], 
  cex = 0.75,
  col = rep(2:4,2), 
  lty =  c(rep(1, 3), rep(2, 3)),
  bty = "n"
) # legend end

##### Summary of intermediate outcomes ----
view(
  dfSummary(
    t(colSums(a_out_interm, dims = 1)), # Convert to data frame automatically
    round.digits = 3,
    style = "grid",
    plain.ascii = FALSE,
    graph.magnif = 1.2,
    headings = FALSE,
    display.labels = FALSE,
    escape.pipe = TRUE,
    varnumbers = FALSE,
    labels.col = FALSE
  )
)
