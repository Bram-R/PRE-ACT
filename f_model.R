f_model <- function(params, intermediate = FALSE) {
  #' State-Transition Model
  #'
  #' This function simulates a state-transition model using a Markov framework. It calculates 
  #' the expected costs, quality-adjusted life years (QALYs), and life-years (LYs) for two 
  #' treatment strategies over a specified time horizon.
  #'
  #' @param params A named list containing model parameters, including:
  #' \itemize{
  #'   \item **Transition Probabilities:**
  #'     \itemize{
  #'       \item \code{tp_ef_ef}: Probability of remaining in "Event-Free".
  #'       \item \code{p_event_lrr}, \code{p_event_death}: Probability of event being locoregional recurrence (LRR) or death.
  #'       \item \code{tp_lrr_dm}, \code{tp_lrr_death}: Probability of transitioning from LRR to distant metastasis (DM) or death.
  #'       \item \code{tp_dm_death}: Probability of transitioning from DM to death.
  #'     }
  #'   \item **Health State Utilities:**
  #'     \itemize{
  #'       \item \code{utility_ef}, \code{utility_lrr}, \code{utility_dm}, \code{utility_death}: QALY weights per state.
  #'     }
  #'   \item **Costs:**
  #'     \itemize{
  #'       \item \code{cost_t1}, \code{cost_t2}: Strategy-specific costs.
  #'       \item \code{cost_ef_y6}, \code{cost_lrr}, \code{cost_dm}, \code{cost_death}: Costs per health state.
  #'       \item \code{cost_lrr_event}, \code{cost_dm_event}, \code{cost_death_event}: Event-related costs.
  #'       \item \code{cost_prev_arm_lymphedema}, \code{cost_prev_pain}, \code{cost_prev_fatigue}, \code{cost_prev_breast_atrophy}: Costs of toxicity prevention.
  #'     }
  #'   \item **Discount Rates:**
  #'     \itemize{
  #'       \item \code{discount_costs}, \code{discount_qalys}, \code{discount_lys}: Discount rates for costs, QALYs, and life-years.
  #'     }
  #' }
  #'
  #' @param intermediate Logical. If \code{FALSE} (default), the function returns a named vector of total 
  #' costs, QALYs, and LYs for each treatment strategy. If \code{TRUE}, it returns detailed 
  #' intermediate cycle-level outputs in a matrix format for further analysis.
  #'
  #' @details
  #' This function:
  #' \itemize{
  #'   \item Defines transition matrices for each treatment.
  #'   \item Ensures transition probabilities are consistent (e.g., adjusting for background mortality).
  #'   \item Estimates state transitions.
  #'   \item Estimates **toxicities over time** using `f_interpolate_toxicity()`.
  #'   \item Computes **discounted costs, QALYs, and LYs** for each treatment strategy.
  #'   \item Returns either a named vector with aggregate results or a matrix with detailed cycle-level outputs.
  #' }
  #'
  #' @return 
  #' If \code{intermediate = FALSE}, returns a named vector containing:
  #' \itemize{
  #'   \item \code{Cost_Current_practice}, \code{Cost_New_treatment}
  #'   \item \code{QALY_Current_practice}, \code{QALY_New_treatment}
  #'   \item \code{LY_Current_practice}, \code{LY_New_treatment}
  #' }
  #'
  #' If \code{intermediate = TRUE}, returns a **matrix** containing:
  #' \itemize{
  #'   \item **Cycle-level health state costs, toxicity costs, and event costs** for each treatment.
  #'   \item **Cycle-level QALYs and LYs** for each treatment.
  #'   \item **Markov trace** (state occupancy probabilities) for each treatment.
  #'   \item **Toxicity prevalence over time** for each treatment.
  #' }
  #'
  #' @examples
  #' # Run model with aggregate results
  #' results <- f_model(f_input(n_sim = 1))
  #' print(results)
  #'
  #' # Run model with intermediate results
  #' results_intermediate <- f_model(f_input(n_sim = 1), intermediate = TRUE)
  #' head(results_intermediate)
  #'
  #' @export
  
  # params <- f_input(n_sim = 1) # for validation purposes
  
  #### Transition matrices ----
  # Initialize transition probability matrices and transition dynamics matrices
  a_transition <- a_transition_dynamics <- array(
    data = 0,
    dim = c(n_treatments, n_t, n_states, n_states),
    dimnames = list(v_treatments, 1:n_t, v_states, v_states)
  )
  
  # Ensure transitions to death are > age and gender matched general population mortality 
  tp_ef_death <- pmax(m_gen_pop_mortality[, 2], rep((1 - params$tp_ef_ef) * params$p_event_death, n_t)) 
  tp_ef_ef <- pmin(params$tp_ef_ef, 1 - tp_ef_death) # correction for consistency (else negative transitions might occur with higher mortality)   
  
  tp_lrr_death <- pmax(m_gen_pop_mortality[, 2], params$tp_lrr_death)
  tp_dm_death <- pmax(m_gen_pop_mortality[, 2], 1 - params$tp_dm_dm)
  
  # Transition probabilities for treatment 1
  # From health state 1 "Event free"
  a_transition[1, , v_states[1], v_states[1]] <- tp_ef_ef                                                 # Stay in health state 1 "Event free"
  a_transition[1, , v_states[1], v_states[2]] <- (1 - tp_ef_ef - tp_ef_death) * params$p_event_lrr        # Transition to health state 2 "Locoregional recurrence"
  a_transition[1, , v_states[1], v_states[3]] <- (1 - tp_ef_ef - tp_ef_death) * (1 - params$p_event_lrr)  # Transition to health state 3 "Distant metastasis"
  a_transition[1, , v_states[1], v_states[4]] <- tp_ef_death                                              # Transition to health state 4 "Death"
  
  # From health state 2 "Locoregional recurrence"
  a_transition[1, , v_states[2], v_states[2]] <- 1 - params$tp_lrr_dm - tp_lrr_death                      # Stay in health state 2 "Locoregional recurrence"
  a_transition[1, , v_states[2], v_states[3]] <- params$tp_lrr_dm                                         # Transition to health state 3 "Distant metastasis"
  a_transition[1, , v_states[2], v_states[4]] <- tp_lrr_death                                             # Transition to health state 4 "Death"
  
  # From health state 3 "Distant metastasis"
  a_transition[1, , v_states[3], v_states[3]] <- 1 - tp_dm_death                                          # Stay in health state 3 "Distant metastasis"
  a_transition[1, , v_states[3], v_states[4]] <- tp_dm_death                                              # Transition to health state 4 "Death"
  
  # From health state 4 "Death"
  a_transition[1, , v_states[4], v_states[4]] <- 1                                                        # Heath state 4 "Death" is absorbing
  
  # Transition probabilities for treatment 2
  a_transition[2, , , ] <- a_transition[1, , , ]                                                          # Copy from treatment 1
  
  # checking whether transitions sum to 1 using testthat package
  # for (i_treatment in 1:n_treatments){ # loop through the treatment options
  #   for (t in 1:n_t){ # loop through the number of cycles
  #     testthat::expect_equal(as.vector(rowSums(a_transition[i_treatment, t, 1:n_states,])), rep(1, n_states))
  #   } # close for loop for cycles
  # } # close for loop for treatments
  
  # check min and max transition probabilities
  # min(a_transition)
  # max(a_transition)
  
  #### Markov trace ----
  # Initialize Markov trace
  a_state_trace <- array(
    data = NA,
    dim = c(n_treatments, n_t + 1, n_states),
    dimnames = list(v_treatments, 0:n_t, v_states)
  )
  a_state_trace[, 1, ] <- matrix(c(1, 0, 0, 0), nrow = n_treatments, ncol = n_states, byrow = TRUE) # Starting health state: 1 "Event free"
  
  # State transitions using nested loops 
  for (i_treatment in 1:n_treatments) {
    for (t in 1:n_t) {
      a_state_trace[i_treatment, t + 1, ] <- a_state_trace[i_treatment, t, ] %*% a_transition[i_treatment, t, , ] # estimate Markov trace for next cycle
      a_transition_dynamics[i_treatment, t, , ] <- diag(a_state_trace[i_treatment, t, ]) %*% a_transition[i_treatment, t, , ] # estimate transition dynamics for next cycle
    }
  }
  
  # checking whether cohort sums to 1 (each cycle) using testthat package
  # for (i_treatment in 1:n_treatments){ # loop through the treatment options
  #     testthat::expect_equal(as.vector(rowSums(a_state_trace[i_treatment, , ])), rep(1, n_t + 1))
  # } # close for loop for treatments
  
  # min(a_state_trace)
  # max(a_state_trace)
  # min(a_transition_dynamics)
  # max(a_transition_dynamics)
  
  #### Toxicity over time ----
  a_toxicity <- array(
    data = NA,
    dim = c(n_treatments, n_t + 1, 4),  
    dimnames = list(treatment = v_treatments, month = 0:n_t, outcome = c("arm_lymphedema", "pain", "fatigue", "breast_atrophy"))
  )
  
  a_toxicity[1,,"arm_lymphedema"] <- f_interpolate_toxicity(data.frame(
    month = c(0, 2, 12, 24, 36, 48, 60, n_t),
    toxicity = c(params$p_arm_lymphedema_m0, params$p_arm_lymphedema_m2, params$p_arm_lymphedema_m12, params$p_arm_lymphedema_m24, 
                 params$p_arm_lymphedema_m36, params$p_arm_lymphedema_m48, params$p_arm_lymphedema_m60, params$p_arm_lymphedema_m60)), 
    monthly_cycles = 0:n_t)[,2]
  
  a_toxicity[2,,"arm_lymphedema"] <- (a_toxicity[1,,"arm_lymphedema"] * params$AI_se_arm_lymphedema) ^ (1 / params$hr_prev_arm_lymphedema) +
    a_toxicity[1,,"arm_lymphedema"] * (1 - params$AI_se_arm_lymphedema) 
  
  a_toxicity[1,,"pain"] <- f_interpolate_toxicity(data.frame(
    month = c(0, 2, 12, 24, 36, 48, 60, n_t),
    toxicity = c(params$p_pain_m0, params$p_pain_m2, params$p_pain_m12, params$p_pain_m24, 
                 params$p_pain_m36, params$p_pain_m48, params$p_pain_m60, params$p_pain_m60)), 
    monthly_cycles = 0:n_t)[,2]
  
  a_toxicity[2,,"pain"] <- a_toxicity[1,,"pain"] * params$AI_se_pain * params$rr_prev_pain +
    a_toxicity[1,,"pain"] * (1 - params$AI_se_pain) 
  
  a_toxicity[1,,"fatigue"] <- f_interpolate_toxicity(data.frame(
    month = c(0, 2, 12, 24, 36, 48, 60, n_t),
    toxicity = c(params$p_fatigue_m0, params$p_fatigue_m2, params$p_fatigue_m12, params$p_fatigue_m24, 
                 params$p_fatigue_m36, params$p_fatigue_m48, params$p_fatigue_m60, params$p_fatigue_m60)), 
    monthly_cycles = 0:n_t)[,2]
  
  a_toxicity[2,,"fatigue"] <- a_toxicity[1,,"fatigue"] * params$AI_se_fatigue * params$rr_prev_fatigue +
    a_toxicity[1,,"fatigue"] * (1 - params$AI_se_fatigue) 
  
  a_toxicity[1,,"breast_atrophy"] <- f_interpolate_toxicity(data.frame(
    month = c(0, 2, 12, 24, 36, 48, 60, n_t),
    toxicity = c(params$p_breast_atrophy_m0, params$p_breast_atrophy_m2, params$p_breast_atrophy_m12, params$p_breast_atrophy_m24, 
                 params$p_breast_atrophy_m36, params$p_breast_atrophy_m48, params$p_breast_atrophy_m60, params$p_breast_atrophy_m60)), 
    monthly_cycles = 0:n_t)[,2]
  
  a_toxicity[2,,"breast_atrophy"] <- a_toxicity[1,,"breast_atrophy"] * params$AI_se_breast_atrophy * params$rr_prev_breast_atrophy +
    a_toxicity[1,,"breast_atrophy"] * (1 - params$AI_se_breast_atrophy) 
  
  # Calculate prevalence per toxicity 
  n_prevalence_arm_lymphedema <- max(a_toxicity[1,,"arm_lymphedema"])
  n_prevalence_pain <- max(a_toxicity[1,,"pain"])
  n_prevalence_fatigue <- max(a_toxicity[1,,"fatigue"])
  n_prevalence_breast_atrophy <- max(a_toxicity[1,,"breast_atrophy"])
  
  # Estimate number of patients receiving preventative measures (i.e. proportion of positive tests)
  n_highrisk_arm_lymphedema <- n_prevalence_arm_lymphedema * params$AI_se_arm_lymphedema + # True positives 
    (1 - n_prevalence_arm_lymphedema) * (1 - params$AI_sp_arm_lymphedema)                  # False positives
  n_highrisk_pain <- n_prevalence_pain * params$AI_se_pain +                               # True positives 
    (1 - n_prevalence_pain) * (1 - params$AI_sp_pain)                                      # False positives
  n_highrisk_fatigue <- n_prevalence_fatigue * params$AI_se_fatigue +                      # True positives 
    (1 - n_prevalence_fatigue) * (1 - params$AI_sp_fatigue)                                # False positives
  n_highrisk_breast_atrophy <- n_prevalence_breast_atrophy * params$AI_se_breast_atrophy + # True positives 
    (1 - n_prevalence_breast_atrophy) * (1 - params$AI_sp_breast_atrophy)                  # False positives
  
  #### Outcomes ----
  # calculate discount weight for each cycle
  m_discount <- matrix(c(1 / (1 + params$discount_costs) ^ (0:n_t / 12), 
                         1 / (1 + params$discount_qalys) ^ (0:n_t / 12),
                         1 / (1 + params$discount_lys) ^ (0:n_t / 12)),
                       nrow = n_t + 1, ncol = 3)
  
  # Cost and (dis)utility matrices
  m_cost <- matrix(c(params$cost_ef_y6,          # Cost for health state 1 "Event free"
                     params$cost_lrr,            # Cost for health state 2 "Locoregional recurrence"
                     params$cost_dm,             # Cost for health state 3 "Distant metastasis"
                     params$cost_death),         # Cost for health state 4 "Death"
                   nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  
  m_cost[1:37, 1] <- params$cost_ef_y1_3         # Time dependent EF health state costs
  m_cost[38:61, 1] <- params$cost_ef_y4_5        # Time dependent EF health state costs   
  
  m_utility <- matrix(c(params$utility_ef,       # Utility for health state 1 "Event free"
                        params$utility_lrr,      # Utility for health state 2 "Locoregional recurrence"
                        params$utility_dm,       # Utility for health state 3 "Distant metastasis"
                        params$utility_death),   # Utility for health state 4 "Death"
                      nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  
  m_utility <- apply(m_utility, 2, function(col) pmin(m_gen_pop_utility[,2], col)) # Ensure utility  are =< age and gender matched general population utility
  
  m_tox_cost <- matrix(c(params$cost_arm_lymphedema,                # Monthly costs for arm_lymphedema
                         params$cost_pain,                          # Monthly costs for pain
                         params$cost_fatigue,                       # Monthly costs for fatigue
                         params$cost_breast_atrophy),               # Monthly costs for breast_atrophy
                       nrow = n_t + 1, ncol = 4, byrow = TRUE)
  
  m_tox_disutility <- matrix(c(params$disutility_arm_lymphedema,    # Disutility for arm_lymphedema
                               params$disutility_pain,              # Disutility for pain
                               params$disutility_fatigue,           # Disutility for fatigue
                               params$disutility_breast_atrophy),   # Disutility for breast_atrophy
                             nrow = n_t + 1, ncol = 4, byrow = TRUE)
  
  # Calculate LYs, QALYs and costs per cycle
  a_costs <- a_state_trace[, , ] * 
    rep(m_cost, each = n_treatments) *                        # Multiply by cost matrix
    rep(m_discount[, 1], each = n_treatments)                 # Multiply by discount factor
  a_costs[1, 1, 1] <- a_costs[1, 1, 1] + params$cost_t1       # Add strategy costs
  a_costs[2, 1, 1] <- a_costs[2, 1, 1] + params$cost_t2       # Add strategy costs
  a_costs_tox <- a_toxicity[, , ] *
    rep(m_tox_cost, each = n_treatments) *                    # Multiply by cost matrix
    rep(m_discount[, 1], each = n_treatments) *               # Multiply by discount factor
    array(1 - a_state_trace[, , 4], dim = c(2, n_t + 1, 4))   # Apply survival component
  
  a_qalys <- a_state_trace[, , ] *
    rep(m_utility, each = n_treatments) *                     # Multiply by utlity matrix
    rep(m_discount[, 2], each = n_treatments) *               # Multiply by discount factor
    1/12                                                      # Multiply by time correction (monthly cycles)
  a_qalys_tox <- a_toxicity[, , ] *
    rep(m_tox_disutility, each = n_treatments) *              # Multiply by disutlity matrix
    rep(m_discount[, 2], each = n_treatments) *               # Multiply by discount factor
    array(1 - a_state_trace[, , 4], dim = c(2, n_t + 1, 4)) * # Apply survival component
    1/12                                                      # Multiply by time correction (monthly cycles)
  
  a_lys <- a_state_trace[, , -4] *                            # Remove health state 4 "Death"  
    rep(m_discount[, 3], each = n_treatments) *               # Multiply by discount factor
    1/12                                                      # Multiply by time correction  (monthly cycles) 
  
  # Event related costs
  m_event_costs <- matrix((a_transition_dynamics[, , 1, 2] * params$cost_lrr_event +      # Costs related to developing locoregional recurrence
                             a_transition_dynamics[, , 1, 3] * params$cost_dm_event +     # Costs related to developing distant metastasis
                             a_transition_dynamics[, , 2, 3] * params$cost_dm_event +     # Costs related to developing distant metastasis
                             a_transition_dynamics[, , 1, 4] * params$cost_death_event +  # Costs related to end of life
                             a_transition_dynamics[, , 2, 4] * params$cost_death_event +  # Costs related to end of life
                             a_transition_dynamics[, , 3, 4] * params$cost_death_event) * # Costs related to end of life
                            rep(m_discount[-1, 1], each = n_treatments),                  # Multiply by discount factor
                          nrow = n_treatments, ncol = n_t)
  
  # Toxicity prevention costs and disutility (assumed in 1st cycle, so no discounting)
  n_tox_prev_costs <- n_highrisk_arm_lymphedema * params$cost_prev_arm_lymphedema +    
    n_highrisk_pain * params$cost_prev_pain + 
    n_highrisk_fatigue * params$cost_prev_fatigue + 
    n_highrisk_breast_atrophy * params$cost_prev_breast_atrophy 
  
  n_tox_prev_disutility <- n_highrisk_arm_lymphedema * params$disutility_prev_arm_lymphedema +    
    n_highrisk_pain * params$disutility_prev_pain + 
    n_highrisk_fatigue * params$disutility_prev_fatigue + 
    n_highrisk_breast_atrophy * params$disutility_prev_breast_atrophy 
  
  #### Results ----
  if(intermediate == FALSE) { 
    v_results <- setNames(
      c(rowSums(a_costs) + rowSums(m_event_costs) + rowSums(a_costs_tox) + c(0, n_tox_prev_costs),  
        rowSums(a_qalys) + rowSums(a_qalys_tox) + c(0, n_tox_prev_disutility), 
        rowSums(a_lys)
      ), # end c     
      c(paste0("Cost_", v_treatments),                          
        paste0("QALY_", v_treatments), 
        paste0("LY_", v_treatments)
      ) # end c
    ) # end setNames
    return(v_results)
  } else {
    # a_res_intermediate <- data.matrix( # convert dataframe to matrix (for storing in array)
    #   data.frame( # create dataframe
    #     "t1.Health.state.costs" = a_costs[1,,],
    #     "t1.Toxicity.costs" = a_costs_tox[1,,],
    #     "t1.Event.costs" = c(0, m_event_costs[1,]),
    #     "t2.Health.state.costs" = a_costs[2,,],
    #     "t2.Toxicity.costs" = a_costs_tox[2,,],
    #     "t2.Event.costs" = c(n_tox_prev_costs, m_event_costs[2,]),
    #     "t1.Health.state.qalys" = a_qalys[1,,],
    #     "t1.Toxicity.qalys" = a_qalys_tox[1,,],
    #     "t1.Event.qalys" = rep(0, n_t + 1),
    #     "t2.Health.state.qalys" = a_qalys[2,,],
    #     "t2.Toxicity.qalys" = a_qalys_tox[2,,],
    #     "t2.Event.qalys" = c(n_tox_prev_disutility, rep(0, n_t)),
    #     "t1.Trace" = a_state_trace[1,,],
    #     "t2.Trace" = a_state_trace[2,,],
    #     "t1.Toxicity" = a_toxicity[1,,],
    #     "t2.Toxicity" = a_toxicity[2,,]
    #   ) # close dataframe
    # ) # close data.matrix
    # return(a_res_intermediate)
    # 
    # create empty array for results per outcome, treatment, cycle and health state
    a_res_intermediate <- array( # array with 4 dimensions (outcome, treatment, cycles, v_states)
      data = NA,
      dim = c( # specify dimensions
        3, # number of outcomes (Costs, QALYs, LYs/toxicity occurrence)
        n_treatments, # number of treatments
        n_t + 1, # number of cycles
        n_states + 5 # number of health states + 1 for event related costs and outcomes + 4 for toxicity 
      ), # close dim
      dimnames = list( # name dimensions
        c("Cost", "QALY", "LY"),
        v_treatments, 
        0:n_t,
        c(v_states, "Event_related", "arm_lymphedema", "pain", "fatigue", "breast_atrophy")
      ) # close dimnames
    ) # close array
    
    a_res_intermediate[1, , , 1:n_states] <- a_costs
    a_res_intermediate[1, , 2:(n_t + 1), (n_states + 1)] <- m_event_costs
    a_res_intermediate[1, 1, 1, (n_states + 1)] <- 0                      # no event costs in first cycle (for current practice)
    a_res_intermediate[1, 2, 1, (n_states + 1)] <- n_tox_prev_costs
    a_res_intermediate[1, , , (n_states + 2):dim(a_res_intermediate)[4]] <- a_costs_tox
    
    a_res_intermediate[2, , , 1:n_states] <- a_qalys
    a_res_intermediate[2, 1, , (n_states + 1)] <- 0                       # no event qalys  (for current practice)
    a_res_intermediate[2, 2, , (n_states + 1)] <-  c(n_tox_prev_disutility, rep(0, n_t))
    a_res_intermediate[2, , , (n_states + 2):dim(a_res_intermediate)[4]] <- a_qalys_tox
    
    a_res_intermediate[3, , , 1:(n_states - 1)] <- a_lys
    a_res_intermediate[3, , , n_states] <- 0                              # no LYs in death health state
    a_res_intermediate[3, , , (n_states + 1)] <- 0                        # no LYs related to events
    a_res_intermediate[3, , , (n_states + 2):dim(a_res_intermediate)[4]] <-  a_toxicity
    
    return(a_res_intermediate)
    # Check if there are any NA values
    # any(is.na(a_res_intermediate))
    # Returns TRUE if at least one NA exists
    
    # Check if there are any negative values
    # any(a_res_intermediate[1,,,] < 0)
    # any(a_res_intermediate[2,,, 1:n_states] < 0) 
    # any(a_res_intermediate[2,,, (n_states + 1):dim(a_res_intermediate)[4]] < 0) # these might be negative (disutilies)
    # any(a_res_intermediate[3,,,] < 0)
    # Returns TRUE if at least one negative value exists
    
  } # close if statement
}

