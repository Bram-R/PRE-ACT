f_input <- function(n_sim = 5000, seed = 12345, setting = 1) {
  #' Generate Inputs for Probabilistic Sensitivity Analysis (PSA)
  #'
  #' This function generates a structured data frame containing input parameters for 
  #' probabilistic sensitivity analysis (PSA). It incorporates various distributions 
  #' (beta, gamma, log-normal, and static values) to account for uncertainty in 
  #' model parameters related to transition probabilities, utilities, and costs.
  #'
  #' @param n_sim Integer. Number of simulations (use `1` for deterministic values). Default is `5000`.
  #' @param seed Integer. Random seed for reproducibility. Default is `12345`.
  #' @param setting Integer. Selected setting `1` = UK, `2` = FR, `3` = NL. Default is `1`.
  #' 
  #' @details
  #' The function:
  #' \itemize{
  #'   \item Sets a random seed to ensure reproducibility.
  #'   \item Adjusts inputs based on setting where appropriate.
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
  df_input_common <- data.frame(
    
    #### Country independent parameters ----
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
    AI_se_arm_lymphedema = generate_sp_se_cor(                           # Based on Canto, HypoG and REQUITE
      mean_sens = 0.850, mean_spec = 0.680, sd_sens = 0.005, sd_spec = 0.011,
      rho = (0.0019*(10^-3)) / sqrt((0.005^2) * (0.011^2)),              # Covariance = 0.0019*(10^-3); mail Guido 20-3-2025
      n = n_sim, is_psa = is_psa, seed = seed)[,1],
    AI_sp_arm_lymphedema = generate_sp_se_cor(                           # Based on Canto, HypoG and REQUITE
      mean_sens = 0.850, mean_spec = 0.680, sd_sens = 0.005, sd_spec = 0.011,
      rho = (0.0019*(10^-3)) / sqrt((0.005^2) * (0.011^2)),              # Covariance = 0.0019*(10^-3); mail Guido 20-3-2025
      n = n_sim, is_psa = is_psa, seed = seed)[,2],
    AI_se_pain = generate_static(0.000, n_sim, is_psa),                  # Treat all as negative
    AI_sp_pain = generate_static(1.000, n_sim, is_psa),                  # Treat all as negative
    AI_se_fatigue = generate_static(0.000, n_sim, is_psa),               # Treat all as negative
    AI_sp_fatigue = generate_static(1.000, n_sim, is_psa),               # Treat all as negative
    AI_se_fibrosis_induration = generate_static(0.000, n_sim, is_psa),        # Treat all as negative
    AI_sp_fibrosis_induration = generate_static(1.000, n_sim, is_psa),        # Treat all as negative    
    
    # Toxicity prevention strategy performance (e.g. related to arm sleeve)
    hr_prev_arm_lymphedema = generate_lognormal(mean = 0.61, ci_low = 0.43, ci_high = 0.85, n_sim, is_psa), # Paramanandam 2022 https://doi.org/10.1200/JCO.21.02567
    rr_prev_pain = generate_static(0, n_sim, is_psa), 
    rr_prev_fatigue = generate_static(0, n_sim, is_psa),                                                    
    rr_prev_fibrosis_induration = generate_static(0, n_sim, is_psa),                                                    
    
    # Toxicity probabilities
    p_arm_lymphedema_m0 = generate_beta(, , n_sim, is_psa, alpha = 45, beta = 1975),                        # Based on REQUITE
    p_arm_lymphedema_m2 = generate_beta(, , n_sim, is_psa, alpha = 58, beta = 1932),                        # Based on REQUITE
    p_arm_lymphedema_m12 = generate_beta(, , n_sim, is_psa, alpha = 60, beta = 1751),                       # Based on REQUITE
    p_arm_lymphedema_m24 = generate_beta(, , n_sim, is_psa, alpha = 66, beta = 1660),                       # Based on REQUITE
    p_arm_lymphedema_m36 = generate_beta(, , n_sim, is_psa, alpha = 21, beta = 579),                        # Based on REQUITE
    p_arm_lymphedema_m48 = generate_beta(, , n_sim, is_psa, alpha = 10, beta = 274),                        # Based on REQUITE
    p_arm_lymphedema_m60 = generate_beta(, , n_sim, is_psa, alpha = 16, beta = 472),                        # Based on REQUITE
    p_arm_lymphedema_m72 = generate_beta(, , n_sim, is_psa, alpha = 2, beta = 43),                          # Based on REQUITE
    
    p_pain_m0  = generate_static(0, n_sim, is_psa),                                                         ## TBD ##   
    p_pain_m2  = generate_static(0, n_sim, is_psa),                                                         ## TBD ##   
    p_pain_m12  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    p_pain_m24  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    p_pain_m36  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    p_pain_m48  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    p_pain_m60  = generate_static(0, n_sim, is_psa),                                                        ## TBD ##   
    p_pain_m72  = generate_static(0, n_sim, is_psa),   
    
    p_fatigue_m0 = generate_static(0, n_sim, is_psa),                                                       ## TBD ## 
    p_fatigue_m2 = generate_static(0, n_sim, is_psa),                                                       ## TBD ## 
    p_fatigue_m12 = generate_static(0, n_sim, is_psa),                                                      ## TBD ## 
    p_fatigue_m24 = generate_static(0, n_sim, is_psa),                                                      ## TBD ##   
    p_fatigue_m36 = generate_static(0, n_sim, is_psa),                                                      ## TBD ##   
    p_fatigue_m48 = generate_static(0, n_sim, is_psa),                                                      ## TBD ##   
    p_fatigue_m60 = generate_static(0, n_sim, is_psa),                                                      ## TBD ##   
    p_fatigue_m72 = generate_static(0, n_sim, is_psa),                                                      ## TBD ##   
    
    p_fibrosis_induration_m0 = generate_static(0, n_sim, is_psa),                                                ## TBD ## 
    p_fibrosis_induration_m2 = generate_static(0, n_sim, is_psa),                                                ## TBD ## 
    p_fibrosis_induration_m12 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    p_fibrosis_induration_m24 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    p_fibrosis_induration_m36 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    p_fibrosis_induration_m48 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    p_fibrosis_induration_m60 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    p_fibrosis_induration_m72 = generate_static(0, n_sim, is_psa),                                               ## TBD ## 
    
    # Health state utilities (Source TA886)
    utility_ef = generate_beta(0.732, 0.021, n_sim, is_psa),                                                # EAG preference in TA886 (Source: Verrill et al 2020)
    utility_lrr = generate_beta((0.732 + 0.603)/2, sqrt(0.021^2 + 0.030^2), n_sim, is_psa),                 # EAG preference in TA886 (assume average of EF and DM utility)
    utility_dm = generate_beta(0.603, 0.030, n_sim, is_psa),                                                # EAG preference in TA886 (Source: Verrill et al 2020)
    utility_death = generate_static(0, n_sim, is_psa),                                                      # Assumption
    
    # Toxicity prevention disutilities (e.g. related to arm sleeve) - one-off
    # disutility_prev_arm_lymphedema = -generate_beta(0.03248, 0.01171, n_sim, is_psa) * 5/12,                #  Table 3 of https://link.springer.com/article/10.1007/s00520-020-05890-3 converted from year to 5 months (arm sleeve for 8h/day from the first day of RT until 3 months after completion of adjuvant radiotherapy = ~ 5 months)      
    disutility_prev_arm_lymphedema = -generate_static(0, n_sim, is_psa),                
    disutility_prev_pain = -generate_static(0, n_sim, is_psa),    
    disutility_prev_fatigue = -generate_static(0, n_sim, is_psa),     
    disutility_prev_fibrosis_induration = -generate_static(0, n_sim, is_psa)                     
  )
  
  #### UK specific parameters ----
  if(setting == 1) { 
    df_input_country_specific <- data.frame(
      # Model settings
      discount_costs = generate_static(0.035, n_sim, is_psa),       # Discount rate for costs (Source: Standard)
      discount_qalys = generate_static(0.035, n_sim, is_psa),       # Discount rate for QALYs (Source: Standard)
      discount_lys = generate_static(0.000, n_sim, is_psa),         # Discount rate for life-years (Source: Standard)
      
      discount_costs_30y = generate_static(0.035, n_sim, is_psa),   # Discount rate for costs (Source: Standard)
      discount_qalys_30y = generate_static(0.035, n_sim, is_psa),   # Discount rate for QALYs (Source: Standard)
      discount_lys_30y = generate_static(0.000, n_sim, is_psa),     # Discount rate for life-years (Source: Standard)
      
      # Toxicity disutilities
      disutility_arm_lymphedema = -generate_beta(0.00948253, 0.0205023, n_sim, is_psa),          # Hypo-G01 analyses            
      disutility_pain = -generate_beta(0.05848352, 0.0193467, n_sim, is_psa),                    # Hypo-G01 analyses
      disutility_fatigue = -generate_beta(0.04044519, 0.0173747, n_sim, is_psa),                 # Hypo-G01 analyses
      disutility_fibrosis_induration = -generate_static(0, n_sim, is_psa),                       # Hypo-G01 analyses
      
      # Strategy costs
      cost_t1 = generate_gamma(2426.842691, 546.4818983, n_sim, is_psa),                           # Costing analysis by Teresa (assumed without WGS)
      cost_t2 = generate_gamma(22.05370746, 5.514792558, n_sim, is_psa),                           # Costing analysis by Teresa (assumed without WGS)  
      
      # Health state costs
      cost_ef_y1_3 = generate_gamma(39.82492114, 39.82492114 * 0.30, n_sim, is_psa),         
      cost_ef_y4_5 = generate_gamma(20.51906533, 20.51906533 * 0.30, n_sim, is_psa),         
      cost_ef_y6 = generate_gamma(3.956117994, 3.956117994 * 0.30, n_sim, is_psa),           
      cost_lrr = generate_gamma(76.48494789, 76.48494789  * 0.30, n_sim, is_psa),          
      cost_dm = generate_gamma(363.9628555, 363.9628555 * 0.30, n_sim, is_psa),            
      cost_death = generate_static(0, n_sim, is_psa),                                              # Assumption
      
      # Progression and mortality costs - one-off 
      cost_lrr_event = generate_gamma(575.9319284, 575.9319284 * 0.30, n_sim, is_psa),             # Used UK HCRU, matched with French unit costs.                           
      cost_dm_event = generate_gamma(575.9319284, 575.9319284 * 0.30, n_sim, is_psa),              # Used UK HCRU, matched with French unit costs.                          
      cost_death_event = generate_gamma(11549.86, 11549.86 * 0.30, n_sim, is_psa),                                    
      
      # Toxicity prevention costs (e.g. arm sleeve) - one-off
      cost_prev_arm_lymphedema_event = generate_gamma(666.3204523, 11.27235178, n_sim, is_psa),    # Used Belgian HCRU and cost, inflated with UK CPI.
      cost_prev_pain_event = generate_static(0, n_sim, is_psa),    
      cost_prev_fatigue_event = generate_static(0, n_sim, is_psa),                                                        
      cost_prev_fibrosis_induration_event = generate_static(0, n_sim, is_psa),                                            
      
      # Toxicity costs (both ongoing and one-off)
      cost_arm_lymphedema = generate_static(0, n_sim, is_psa),        
      cost_pain = generate_gamma(26.07, 26.07 * 0.3, n_sim, is_psa),                                                             
      cost_fatigue = generate_static(0, n_sim, is_psa),         
      cost_fibrosis_induration = generate_gamma(26.07, 26.07 * 0.3, n_sim, is_psa),                                                  
      
      cost_arm_lymphedema_event = generate_gamma(2448.8685, 2448.8685 * 0.3, n_sim, is_psa),  
      cost_pain_event = generate_gamma(1504.09, 1504.09 * 0.3, n_sim, is_psa),  
      cost_fatigue_event = generate_gamma(2325.55, 2325.55 * 0.3, n_sim, is_psa),           
      cost_fibrosis_induration_event = generate_gamma(34.59, 34.59 * 0.3, n_sim, is_psa)                                  
    )
  }
  
  #### FR specific parameters ----
  if(setting == 2) { 
    df_input_country_specific <- data.frame(
      # Model settings
      discount_costs = generate_static(0.025, n_sim, is_psa),       # Discount rate for costs (Source: Standard)
      discount_qalys = generate_static(0.025, n_sim, is_psa),       # Discount rate for QALYs (Source: Standard)
      discount_lys = generate_static(0.000, n_sim, is_psa),         # Discount rate for life-years (Source: Standard)
      
      discount_costs_30y = generate_static(0.015, n_sim, is_psa),   # Discount rate for costs (Source: Standard)
      discount_qalys_30y = generate_static(0.015, n_sim, is_psa),   # Discount rate for QALYs (Source: Standard)
      discount_lys_30y = generate_static(0.000, n_sim, is_psa),     # Discount rate for life-years (Source: Standard)
      
      # Toxicity disutilities
      disutility_arm_lymphedema = -generate_beta(0.0182440, 0.0229990, n_sim, is_psa),          # Hypo-G01 analyses            
      disutility_pain = -generate_beta(0.0513998, 0.0231084, n_sim, is_psa),                    # Hypo-G01 analyses
      disutility_fatigue = -generate_beta(0.0376353, 0.023412, n_sim, is_psa),                  # Hypo-G01 analyses
      disutility_fibrosis_induration = -generate_static(0, n_sim, is_psa),                      # Hypo-G01 analyses

      # Strategy costs
      cost_t1 = generate_gamma(5510.678541, 2486.447128, n_sim, is_psa),                           # Costing analysis by Teresa (assumed without WGS)
      cost_t2 = generate_gamma(36.2043966, 36.80058921, n_sim, is_psa),                            # Costing analysis by Teresa (assumed without WGS)  
      
      # Health state costs
      cost_ef_y1_3 = generate_gamma(16.67481656, 16.67481656 * 0.30, n_sim, is_psa),         
      cost_ef_y4_5 = generate_gamma(10.54076903, 10.54076903 * 0.30, n_sim, is_psa),         
      cost_ef_y6 = generate_gamma(2.93526415, 2.93526415 * 0.30, n_sim, is_psa),           
      cost_lrr = generate_gamma(20.96592963, 20.96592963 * 0.30, n_sim, is_psa),          
      cost_dm = generate_gamma(70.68341114, 70.68341114 * 0.30, n_sim, is_psa),            
      cost_death = generate_static(0, n_sim, is_psa),                                              # Assumption
      
      # Progression and mortality costs - one-off 
      cost_lrr_event = generate_gamma(140.6070584, 140.6070584 * 0.30, n_sim, is_psa),             # Used UK HCRU, matched with French unit costs.                           
      cost_dm_event = generate_gamma(140.6070584, 140.6070584 * 0.30, n_sim, is_psa),              # Used UK HCRU, matched with French unit costs.                          
      cost_death_event = generate_gamma(4803.938484, 4803.938484 * 0.30, n_sim, is_psa),                                    
      
      # Toxicity prevention costs (e.g. arm sleeve) - one-off
      cost_prev_arm_lymphedema_event = generate_gamma(652.4742443, 11.03811114, n_sim, is_psa),    # Used Belgian HCRU and cost, inflated with French CPI.
      cost_prev_pain_event = generate_static(0, n_sim, is_psa),                                                           
      cost_prev_fatigue_event = generate_static(0, n_sim, is_psa),                                                        
      cost_prev_fibrosis_induration_event = generate_static(0, n_sim, is_psa),                                            
      
      # Toxicity costs (both ongoing and one-off)
      cost_arm_lymphedema = generate_static(0, n_sim, is_psa),        
      cost_pain = generate_static(0, n_sim, is_psa),                                                                
      cost_fatigue = generate_static(0, n_sim, is_psa),               
      cost_fibrosis_induration = generate_gamma(11.11443596, 11.11443596 * 0.3, n_sim, is_psa),                                                  
      
      cost_arm_lymphedema_event = generate_gamma(1358.11, 1358.11 * 0.3, n_sim, is_psa),  
      cost_pain_event = generate_gamma(1565.89, 1565.89 * 0.3, n_sim, is_psa),                                                          
      cost_fatigue_event = generate_gamma(1461.09, 1461.09 * 0.3, n_sim, is_psa),    
      cost_fibrosis_induration_event = generate_gamma(40.53, 40.53 * 0.3, n_sim, is_psa)                                           
    )
  }
  
  #### NL specific parameters ----
  if(setting == 3) { 
    df_input_country_specific <- data.frame(
      # Model settings
      discount_costs = generate_static(0.030, n_sim, is_psa),       # Discount rate for costs (Source: Standard)
      discount_qalys = generate_static(0.015, n_sim, is_psa),       # Discount rate for QALYs (Source: Standard)
      discount_lys = generate_static(0.000, n_sim, is_psa),         # Discount rate for life-years (Source: Standard)
      
      discount_costs_30y = generate_static(0.030, n_sim, is_psa),   # Discount rate for costs (Source: Standard)
      discount_qalys_30y = generate_static(0.015, n_sim, is_psa),   # Discount rate for QALYs (Source: Standard)
      discount_lys_30y = generate_static(0.000, n_sim, is_psa),     # Discount rate for life-years (Source: Standard)
      
      # Toxicity disutilities
      disutility_arm_lymphedema = -generate_beta(0.00268582, 0.0240935, n_sim, is_psa),          # Hypo-G01 analyses            
      disutility_pain = -generate_beta(0.05164852, 0.0193217, n_sim, is_psa),                    # Hypo-G01 analyses
      disutility_fatigue = -generate_beta(0.03450248, 0.0156650, n_sim, is_psa),                 # Hypo-G01 analyses
      disutility_fibrosis_induration = -generate_static(0, n_sim, is_psa),                       # Hypo-G01 analyses
      
      # Strategy costs
      cost_t1 = generate_gamma(18545.1022, 4451.567957, n_sim, is_psa),                            # Costing analysis by Teresa (assumed without WGS)
      cost_t2 = generate_gamma(37.65182529, 25.08999474, n_sim, is_psa),                           # Costing analysis by Teresa (assumed without WGS)  
      
      # Health state costs
      cost_ef_y1_3 = generate_gamma(37.73262451, 37.73262451 * 0.30, n_sim, is_psa),         
      cost_ef_y4_5 = generate_gamma(22.65858917, 22.65858917 * 0.30, n_sim, is_psa),         
      cost_ef_y6 = generate_gamma(3.99765153, 3.99765153 * 0.30, n_sim, is_psa),           
      cost_lrr = generate_gamma(90.56512927, 90.56512927  * 0.30, n_sim, is_psa),          
      cost_dm = generate_gamma(301.90230139, 301.90230139 * 0.30, n_sim, is_psa),            
      cost_death = generate_static(0, n_sim, is_psa),                                              # Assumption
      
      # Progression and mortality costs - one-off 
      cost_lrr_event = generate_gamma(660.1161946, 660.1161946 * 0.30, n_sim, is_psa),             # Used UK HCRU, matched with Dutch unit costs.                           
      cost_dm_event = generate_gamma(660.1161946, 660.1161946 * 0.30, n_sim, is_psa),              # Used UK HCRU, matched with Dutch unit costs.                          
      cost_death_event = generate_gamma(6183.78592, 6183.78592 * 0.30, n_sim, is_psa),                                    
      
      # Toxicity prevention costs (e.g. arm sleeve) - one-off
      cost_prev_arm_lymphedema_event = generate_gamma(805.311062, 13.62369945, n_sim, is_psa),    # Used Belgian HCRU and cost, inflated with Dutch CPI.
      cost_prev_pain_event = generate_static(0, n_sim, is_psa),                                                           
      cost_prev_fatigue_event = generate_static(0, n_sim, is_psa),                                                        
      cost_prev_fibrosis_induration_event = generate_static(0, n_sim, is_psa),                                            
      
      # Toxicity costs (both ongoing and one-off)
      cost_arm_lymphedema = generate_static(0, n_sim, is_psa),        
      cost_pain = generate_gamma(22.51, 22.51 * 0.3, n_sim, is_psa),                                                                    
      cost_fatigue = generate_gamma(26.85, 26.85 * 0.3, n_sim, is_psa),                
      cost_fibrosis_induration = generate_gamma(15.77411819, 15.77411819 * 0.3, n_sim, is_psa),                                                  
      
      cost_arm_lymphedema_event = generate_gamma(1927.49016, 1927.49016 * 0.3, n_sim, is_psa),  
      cost_pain_event = generate_gamma(794.40, 794.40 * 0.3, n_sim, is_psa), 
      cost_fatigue_event = generate_gamma(172.14, 172.14 * 0.3, n_sim, is_psa),           
      cost_fibrosis_induration_event = generate_gamma(40.53, 40.53 * 0.3, n_sim, is_psa)                                           
    )
  }
  
  # Safeguard for invalid setting value
  if (!exists("df_input_country_specific")) {
    stop("Invalid setting: must be 1 (UK), 2 (FR), or 3 (NL)")
  }
  
  # Assembly of df_input
  df_input <- data.frame(
    df_input_common,
    df_input_country_specific
  )
  
  # Dependent calculations
  df_input$tp_lrr_death <- (1 - df_input$tp_ef_ef) * df_input$p_event_death
  df_input$cost_t2 <- df_input$cost_t1 + df_input$cost_t2
  
  # df_input$disutility_arm_lymphedema <- pmin(df_input$utility_arm_lymphedema - df_input$utility_ef, 0)
  # df_input <- subset(df_input, select = -utility_arm_lymphedema)                                                  # remove temporary variable utility_arm_lymphedema
  
  return(df_input)
}
