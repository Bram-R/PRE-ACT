f_input <- function(n_sim = 10000, seed = 12345) {
  #' Generate Inputs for Probabilistic Sensitivity Analysis (PSA)
  #'
  #' This function generates a structured data frame containing input parameters for 
  #' probabilistic sensitivity analysis (PSA). It incorporates various distributions 
  #' (beta, gamma, log-normal, and static values) to account for uncertainty in 
  #' model parameters related to transition probabilities, utilities, and costs.
  #'
  #' @param n_sim Integer. Number of simulations (use `1` for deterministic values). Default is `10000`.
  #' @param seed Integer. Random seed for reproducibility. Default is `12345`.
  #' 
  #' @details
  #' The function:
  #' \itemize{
  #'   \item Sets a random seed to ensure reproducibility.
  #'   \item Uses helper functions (`generate_static`, `generate_beta`, `generate_gamma`, `generate_lognormal`)
  #'         to sample values from different probability distributions.
  #'   \item Defines transition probabilities, health state utilities, toxicity parameters, and costs.
  #'   \item Adjusts probability dependencies where needed (e.g., `tp_lrr_death` based on other inputs).
  #' }
  #' 
  #' The generated inputs align with the structure needed for state-transition modeling in economic evaluations.
  #'
  #' @return A data frame (`df_input`) where each row corresponds to a unique set of sampled PSA parameters.
  #'
  #' @examples
  #' # Generate a small set of PSA inputs (100 simulations)
  #' psa_inputs <- f_input(n_sim = 100, seed = 12345)
  #'
  #' # Generate deterministic inputs (single set of parameters)
  #' deterministic_inputs <- f_input(n_sim = 1, seed = 12345)
  #'
  #' @export
  
  source("f_input_helper.R") # Load helper functions
  
  set.seed(seed)
  is_psa <- n_sim > 1
  
  # Parameter data frame
  df_input <- data.frame(
    
    # Model settings
    discount_costs = generate_static(0.035, n_sim, is_psa),       # Discount rate for costs (Source: Standard)
    discount_qalys = generate_static(0.035, n_sim, is_psa),       # Discount rate for QALYs (Source: Standard)
    discount_lys = generate_static(0.000, n_sim, is_psa),         # Discount rate for life-years (Source: Standard)
    
    # State-transition model probabilities
    # from event free (Source: KEYNOTE-522, NICE TA851)
    tp_ef_ef = generate_beta(0.835^(1/42), (0.805^(1/42) - 0.860^(1/42)) / (2 * qnorm(.975)), n_sim, is_psa),  # Probability to stay in EF; 1 - Event free survival (42m EFS: 83.5% in PCP arm (Source: KEYNOTE-522, NICE TA851)
    p_event_lrr = generate_beta(, , n_sim, is_psa, alpha = 38 + 31, beta = 85 + 62),                           # Probability of EFS event being locoregional recurrence based on both arms (Source: KEYNOTE-522, NICE TA851)
    p_event_death = generate_beta(, , n_sim, is_psa, alpha = 15 + 6, beta = 108 + 87),                         # Probability of EFS event being death based on both arms (Source: KEYNOTE-522, NICE TA851)
    
    # from locoregional recurrence (Source: Hamilton SN, Tyldesley S, Li D, et al. Second malignancies after adjuvant radiation therapy for early stage breast cancer: is there increased risk with addition of regional radiation to local radiation? Int J Radiat Oncol Biol Phys 2015;91:977-985)
    tp_lrr_dm = generate_beta(1 - (1 - 0.73)^(1/60), (1 - (1 - 0.73)^(1/60)) * 0.30, n_sim, is_psa),           # Probability to transit from lrr to dm (Source: Nielsen 2006 doi:10.1016/j.radonc.2006.04.006 Table 4)
    # tp_lrr_dm = generate_beta(0.00757, 0.000012, n_sim, is_psa),                                             # Probability to transit from lrr to dm (Source: clarification response Table 17 of TA424)
    tp_lrr_death = generate_static(NA, n_sim, is_psa),                                                         # Probability to transit from lrr to death (assume the same as from ef to death, calculated below). This assumption is based on TA851 stating that the proportion of LRR -> death is very small and has minimal impact 
    
    # from distant metastasis (Source: KEYNOTE-355, +/- 95% had metastatic disease and +/- locoregional disease)
    tp_dm_dm = generate_beta(0.478^(1/18), (0.436^(1/18) - 0.519^(1/18)) / (2 * qnorm(.975)), n_sim, is_psa),  # OS (18M OS: 47.8% in PCP arm (Source: KEYNOTE-522, NICE TA851)
    
    # Diagnostic performance AI
    AI_se_arm_lymphedema = generate_beta(0.850, 0.005, n_sim, is_psa),   # Based on Canto, HypoG and REQUITE
    AI_sp_arm_lymphedema = generate_beta(0.680, 0.011, n_sim, is_psa),   # Based on Canto, HypoG and REQUITE
    AI_se_breast_atrophy = generate_static(0.000, n_sim, is_psa),        # Treat all as negative
    AI_sp_breast_atrophy = generate_static(1.000, n_sim, is_psa),        # Treat all as negative    
    AI_se_fatigue = generate_static(0.000, n_sim, is_psa),               # Treat all as negative
    AI_sp_fatigue = generate_static(1.000, n_sim, is_psa),               # Treat all as negative
    AI_se_pain = generate_static(0.000, n_sim, is_psa),                  # Treat all as negative
    AI_sp_pain = generate_static(1.000, n_sim, is_psa),                  # Treat all as negative
    
    # Toxicity prevention strategy performance (e.g. related to arm sleeve)
    hr_prev_arm_lymphedema = generate_lognormal(mean = 0.61, ci_low = 0.43, ci_high = 0.85, n_sim, is_psa), # Paramanandam 2022 https://doi.org/10.1200/JCO.21.02567
    rr_prev_breast_atrophy = generate_static(0, n_sim, is_psa),                                             ## TBD ## 
    rr_prev_fatigue = generate_static(0, n_sim, is_psa),                                                    ## TBD ## 
    rr_prev_pain = generate_static(0, n_sim, is_psa),                                                       ## TBD ## 
    
    # Toxicity probabilities
    p_arm_lymphedema_m0 = generate_beta(, , n_sim, is_psa, alpha = 49 + 45, beta = 1205 + 1975),            # Based on HypoG and REQUITE
    p_arm_lymphedema_m2 = generate_beta(, , n_sim, is_psa, alpha = 89 + 58, beta = 1104 + 1932),            # Based on HypoG and REQUITE
    p_arm_lymphedema_m12 = generate_beta(, , n_sim, is_psa, alpha = 130 + 211, beta = 1070 + 1751),         # Based on HypoG and REQUITE
    p_arm_lymphedema_m24 = generate_beta(, , n_sim, is_psa, alpha = 120 + 296, beta = 1066 + 1660),         # Based on HypoG and REQUITE
    p_arm_lymphedema_m36 = generate_beta(, , n_sim, is_psa, alpha = 113 + 579, beta = 1082 + 1422),         # Based on HypoG and REQUITE
    p_arm_lymphedema_m48 = generate_beta(, , n_sim, is_psa, alpha = 60 + 274, beta = 1165 + 1738),          # Based on HypoG and REQUITE
    p_arm_lymphedema_m60 = generate_beta(, , n_sim, is_psa, alpha = 23 + 472, beta = 1229 + 1534),          # Based on HypoG and REQUITE
    
    p_breast_atrophy_m0 = generate_static(0, n_sim, is_psa),                                                ## TBD ## 
    p_breast_atrophy_m2 = generate_static(0, n_sim, is_psa),                                                ## TBD ## 
    p_breast_atrophy_m12 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    p_breast_atrophy_m24 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    p_breast_atrophy_m36 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    p_breast_atrophy_m48 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    p_breast_atrophy_m60 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    
    p_fatigue_m0 = generate_static(0, n_sim, is_psa),                                                       ## TBD ## 
    p_fatigue_m2 = generate_static(0, n_sim, is_psa),                                                       ## TBD ## 
    p_fatigue_m12 = generate_static(0, n_sim, is_psa),                                                      ## TBD ## 
    p_fatigue_m24 = generate_static(0, n_sim, is_psa),                                                      ## TBD ##   
    p_fatigue_m36 = generate_static(0, n_sim, is_psa),                                                      ## TBD ##   
    p_fatigue_m48 = generate_static(0, n_sim, is_psa),                                                      ## TBD ##   
    p_fatigue_m60 = generate_static(0, n_sim, is_psa),                                                      ## TBD ##   
    
    p_pain_m0  = generate_static(0, n_sim, is_psa),                                                         ## TBD ##   
    p_pain_m2  = generate_static(0, n_sim, is_psa),                                                         ## TBD ##   
    p_pain_m12  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    p_pain_m24  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    p_pain_m36  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    p_pain_m48  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    p_pain_m60  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    
    # Health state utilities (Source TA886)
    utility_ef = generate_beta(0.732, 0.021, n_sim, is_psa),                                     # EAG preference in TA886 (Source: Verrill et al 2020)
    utility_lrr = generate_beta((0.732 + 0.603)/2, sqrt(0.021^2 + 0.030^2), n_sim, is_psa),      # EAG preference in TA886 (assume average of EF and DM utility)
    utility_dm = generate_beta(0.603, 0.030, n_sim, is_psa),                                     # EAG preference in TA886 (Source: Verrill et al 2020)
    utility_death = generate_static(0, n_sim, is_psa),                                           # Assumption
    
    # Toxicity prevention disutilities (e.g. related to arm sleeve) - one-off
    disutility_prev_arm_lymphedema = -generate_beta(0.03248, 0.01171, n_sim, is_psa)* 5/12,      # Assumed equal to disutility for rash (Nafees et al 2008 https://pmc.ncbi.nlm.nih.gov/articles/PMC2579282/) for a duration of 5 months (arm sleeve for 8h/day from the first day of RT until 3 months after completion of adjuvant radiotherapy = ~ 5 months)
    disutility_prev_breast_atrophy = -generate_static(0, n_sim, is_psa),                         ## TBD ##     
    disutility_prev_fatigue = -generate_static(0, n_sim, is_psa),                                ## TBD ## 
    disutility_prev_pain = -generate_static(0, n_sim, is_psa),                                   ## TBD ## 
    
    # Toxicity disutilities
    utility_arm_lymphedema = generate_beta(0.633, 0.197, n_sim, is_psa),                         # Utility reported for lymphedema patients in https://onlinelibrary.wiley.com/doi/full/10.1111/wrr.12529  
    disutility_arm_lymphedema = generate_static(NA, n_sim, is_psa),                              # Assumed mean EF utility - 0.633 (calculated below)  
    disutility_breast_atrophy = -generate_static(0, n_sim, is_psa),                              ## TBD ##  
    disutility_fatigue = -generate_beta(0.115, 0.115, n_sim, is_psa),                            # TA612 CS Table 53
    disutility_pain = -generate_static(0, n_sim, is_psa),                                        ## TBD ## 
    
    # Health state costs (Source TA851) --> Year: 2020
    cost_ef_y1_3 = generate_gamma(7.55, 7.55 * 0.30, n_sim, is_psa) * (365.25 / 7 / 12),         # TA851 CS Table 62 converted from week to month (Source: TA424)
    cost_ef_y4_5 = generate_gamma(3.89, 3.89 * 0.30, n_sim, is_psa) * (365.25 / 7 / 12),         # TA851 CS Table 62 converted from week to month (Source: TA424)
    cost_ef_y6 = generate_gamma(0.75, 0.75 * 0.30, n_sim, is_psa) * (365.25 / 7 / 12),           # TA851 CS Table 62 converted from week to month (Source: TA424) 
    cost_lrr = generate_gamma(14.50, 14.50  * 0.30, n_sim, is_psa) * (365.25 / 7 / 12),          # TA851 CS Table 62 converted from week to month (Source: TA569)
    cost_dm = generate_gamma(69.00, 69.00 * 0.30, n_sim, is_psa) * (365.25 / 7 / 12),            # TA851 CS Table 62 converted from week to month (Source: TA546)
    cost_death = generate_static(0, n_sim, is_psa),                                              # Assumption
    
    # Strategy costs
    cost_t1 = generate_gamma(2590.22, 641.97, n_sim, is_psa),                                  # Costing analysis by Teresa (assumed without WGS)
    cost_t2 = generate_gamma(2612.74, 642.09, n_sim, is_psa),                                  # Costing analysis by Teresa (assumed without WGS)  
    
    # Toxicity prevention costs (e.g. arm sleeve) - one-off --> Year: 2020
    cost_prev_arm_lymphedema = generate_gamma(521.29 + 122.29, sqrt(105.83^2 + 105.85^2), n_sim, is_psa) * 5/12,  #  Table 3 of https://link.springer.com/article/10.1007/s00520-020-05890-3 converted from year to 5 months (arm sleeve for 8h/day from the first day of RT until 3 months after completion of adjuvant radiotherapy = ~ 5 months)
    cost_prev_breast_atrophy = generate_static(0, n_sim, is_psa),                                                 ## TBD ##    
    cost_prev_fatigue = generate_static(0, n_sim, is_psa),                                                        ## TBD ##    
    cost_prev_pain = generate_static(0, n_sim, is_psa),                                                           ## TBD ## 
    
    # Event costs (Source TA851) --> Year: 2020 
    cost_lrr_event = generate_gamma(474.76, 474.76 * 0.30, n_sim, is_psa),                                        # costs of oncologist visit, CT scan, full blood count and MRI scan
    cost_dm_event = generate_gamma(474.76, 474.76 * 0.30, n_sim, is_psa),                                         # costs of oncologist visit, CT scan, full blood count and MRI scan 
    cost_death_event = generate_gamma(8347.03, 8347.03 * 0.30, n_sim, is_psa),                                    # costs of district nurse, nursing & residential care, hospice care and Marie Curie nursing service (Georgiou et al 2014)
    
    # Toxicity costs --> Year: 2020 for cost_arm_lymphedema and 2018 for cost_fatigue
    cost_arm_lymphedema = generate_gamma(1803.35 + 445.58, sqrt(343.48^2 + 246.40^2), n_sim, is_psa) / 12,        # Table 3 of https://link.springer.com/article/10.1007/s00520-020-05890-3
    cost_breast_atrophy = generate_static(0, n_sim, is_psa),                                                      ## TBD ## 
    cost_fatigue = generate_gamma(3000.39, 300.03, n_sim, is_psa) / (9.43 * 7) * (365.25 / 12),                   # TA612 CS Table 53 (2018 price year), duration: 9.43 weeks CS Table 40
    cost_pain = generate_static(0, n_sim, is_psa)                                                                 ## TBD ## 
  )
  
  df_input$tp_lrr_death <- (1 - df_input$tp_ef_ef) * df_input$p_event_death
  df_input$disutility_arm_lymphedema <- pmin(df_input$utility_arm_lymphedema - df_input$utility_ef, 0)
  df_input <- subset(df_input, select = -utility_arm_lymphedema)                                                  # remove temporary variable utility_arm_lymphedema
  
  return(df_input)
}
