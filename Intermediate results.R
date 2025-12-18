#### Description: Intermediate results of a State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# Clear workspace
rm(list = ls())

# Load custom functions
source("Model setup.R")            # Model setup and definitions

#### Model Inputs ----
# Country specific
n_currency <- "Pound"
n_wtp <- 30000 # willingness to pay value
v_wtp <- seq(from = 0, to = 50000, by = 2000)  # willingness to pay vector

# Create a dataframe for probabilistic sensitivity analysis (PSA) inputs
m_gen_pop_utility <- f_gen_pop_utility(n_age_baseline = n_age_baseline, n_t = n_t)
m_gen_pop_mortality <- f_gen_pop_mortality(n_age_baseline = n_age_baseline, n_t = n_t, n_p_female = n_p_female)
df_input <- f_input(n_sim = n_sim)

##### Obtain intermediate results using wrapper intermediate function ----
# Intermediate outcomes are stored in a 3D array with the following structure: cycles, outcomes, simulations. 
# Here outcomes are Costs, QALYs and LYs/incidence are further specified by:
# - Health state: Costs, QALYs and LYs associated with health state occupancy.
# - Toxicity: Costs, QALYs and incidence associated with the occurrence of toxicities (either decreased HRQoL, costs related to toxicity management or toxicity incidence).
# - Event: One-off costs and QALYs, i.e. treatment costs, toxicity prevention costs (e.g. arm sleeve) as well as one-off costs related to the development of recurrence (either loco-regional or distant) as well as mortality. 

# Probabilistic results (intermediate)
a_out_interm <- f_wrapper_intermediate(df_input)

# str(a_out_interm[,,]) # inspect
# names(a_out_interm[1, , 1]) # inspect

##### Summary of intermediate outcomes ----
v_res <- rowMeans(colSums(a_out_interm, dims = 1))

v_res

# Validate intermediate outcomes (comparing the results below with the probabilistic base-case results)
# sum(v_res[1:9])
# sum(v_res[10:18])
# sum(v_res[19:27])
# sum(v_res[28:36])

txt <- capture.output(
  print(
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
)
writeLines(txt, "text/Intermediate_results.txt")

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
png(file = paste0("plots/", "costs_vs_time_", v_treatments[1], ".png"), width = 1500, height = 1500)
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
dev.off()

png(file = paste0("plots/", "costs_vs_time_", v_treatments[2], ".png"), width = 1500, height = 1500)
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
dev.off()

png(file = paste0("plots/", "costs_vs_time", ".png"), width = 1500, height = 1500)
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
dev.off()

# Probabilistic QALYs 
png(file = paste0("plots/", "qalys_vs_time_", v_treatments[1], ".png"), width = 1500, height = 1500)
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
dev.off()

png(file = paste0("plots/", "qalys_vs_time_", v_treatments[2], ".png"), width = 1500, height = 1500)
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
dev.off()

png(file = paste0("plots/", "qalys_vs_time", ".png"), width = 1500, height = 1500)
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
dev.off()

##### Toxicity per cycle ----
m_toxicity <- rowMeans(a_out_interm[, c(41:44, 50:53), ], dims = 2)
m_toxicity_percentiles <- cbind(
  rowQuantiles(a_out_interm[, 41, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 42, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 43, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 44, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 50, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 51, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 52, ], probs = c(0.025, 0.975)),
  rowQuantiles(a_out_interm[, 53, ], probs = c(0.025, 0.975))
)

# Toxicity incidence
png(file = paste0("plots/", "tox_vs_time", ".png"), width = 1500, height = 1500)
matplot(
  x = 0:n_t,
  y = m_toxicity[, 1:8], 
  ylim = c(0, 0.1),
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
dev.off()

# Arm lymphedema incidence
png(file = paste0("plots/", "arm_lymphedema_vs_time", ".png"), width = 1500, height = 1500)
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
  c(paste0("Probabilistic mean ", v_treatments[1]), dimnames(m_toxicity_percentiles[, 1:2])[[2]],
    paste0("Probabilistic mean ", v_treatments[2]), dimnames(m_toxicity_percentiles[, 1:2])[[2]]),
  cex = 0.5,
  col = c(2, 2, 2, 3, 3, 3),
  lty = c(1, 2, 2, 1, 2, 2),
  bty = "n"
)
dev.off()

##### LYs cycle ----
m_trace <- rowMeans(a_out_interm[, c(37:40, 46:49), ], dims = 2)

# Probabilistic LYs
png(file = paste0("plots/", "ly_by_hs_vs_time", ".png"), width = 1500, height = 1500)
matplot(
  x = 0:n_t,
  y = m_trace[, 1:8], 
  #ylim = c(0, 1),
  type = "l",
  ylab = "LYs",
  xlab = "Cycle",
  main = paste0("Average LYs for ", v_treatments[1], " and ", v_treatments[2]),
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
dev.off()

##### Disaggregated costs and QALYs per cycle ----
m_out_dis_cycle_res <- rowMeans(a_out_interm[, 1:36, ], dims = 2)

# Health state costs
png(file = paste0("plots/", "costs_by_hs_vs_time", ".png"), width = 1500, height = 1500)
matplot(
  x = 0:n_t,
  y = m_out_dis_cycle_res[, 1:18], 
  type = "l",
  ylab = "Costs",
  xlab = "Cycle",
  main = paste0("Average health state costs for ", v_treatments[1], " and ", v_treatments[2]),
  col = rainbow(9),
  lty = c(rep(1, 9), rep(2, 9))
) # matplot end
legend(
  "topright", 
  inset = c(0.05, 0),
  dimnames(m_out_dis_cycle_res[, 1:18])[[2]], 
  cex = 0.75,
  col = rainbow(9), 
  lty =  c(rep(1, 9), rep(2, 9)),
  bty = "n"
) # legend end
dev.off()

# Health state QALYs
png(file = paste0("plots/", "qalys_by_hs_vs_time", ".png"), width = 1500, height = 1500)
matplot(
  x = 0:n_t,
  y = m_out_dis_cycle_res[, 19:36], 
  #ylim = c(0, 1),
  type = "l",
  ylab = "Costs",
  xlab = "Cycle",
  main = paste0("Average health state QALYs ", v_treatments[1], " and ", v_treatments[2]),
  col = rainbow(9),
  lty = c(rep(1, 9), rep(2, 9))
) # matplot end
legend(
  "topright", 
  inset = c(0.05, 0),
  dimnames(m_out_dis_cycle_res[, 19:36])[[2]], 
  cex = 0.75,
  col = rainbow(9), 
  lty =  c(rep(1, 9), rep(2, 9)),
  bty = "n"
) # legend end
dev.off()
