###################################################################################################
#Title: Model set-up (run script for each new parameter set)
#Author: Hailey Park
#Date: April 18, 2024
###################################################################################################

#Load in data
severe_waning_data <- read.csv("data/clean-data/combined_severe_waning_predictions_weekly.csv")[,-1]
nonsevere_waning_data <- read.csv("data/clean-data/combined_nonsevere_waning_predictions_weekly.csv")[,-1] 
severe_waning_prior_inf_only_data <- read.csv("data/clean-data/combined_severe_waning_predictions_weekly_prior_inf_only.csv")[,-1]
nonsevere_waning_prior_inf_only_data <- read.csv("data/clean-data/combined_nonsevere_waning_predictions_weekly_prior_inf_only.csv")[,-1]
average_severe_incidence <- read.csv("data/clean-data/observed-inc-estimates.csv")[,-1]

model_set_up_fn <- function(params){
  #Set varied params (lambda, pe vars, nonsevere mult, etc.)
  lambda <- params[1]
  pe_nonsevere_level_vacc <- params[2]
  pe_nonsevere_level_inf <- params[3]
  pe_nonsevere_level_hybrid <- params[4]
  pe_nonsevere_rate_vacc <- params[5]
  pe_nonsevere_rate_inf <- params[6]
  pe_nonsevere_rate_hybrid <- params[7]
  nonsevere_mult <- params[8]
  
  
  #Clean waning curve data and add adjustments (if rate waning > 20%, revisit formula)
  severe_waning <- severe_waning_data %>% 
    filter(estimate == 'mean') %>% 
    dplyr::select(-c(study, estimate, month_input, months)) %>%
    rowwise() %>% 
    mutate(ve_pred = max(ve_pred, 0)) %>%
    rename(severe_ve_pred = ve_pred)
  
  severe_waning_prior_inf_only <- severe_waning_prior_inf_only_data %>%
    filter(estimate == "mean") %>%
    rowwise() %>%
    mutate(ve_pred = max(ve_pred, 0)) %>%
    rename(severe_prior_inf_only_ve_pred = ve_pred) %>%
    dplyr::select(age_group, weeks, immunocompromised, severe_prior_inf_only_ve_pred)
  
  nonsevere_waning <- nonsevere_waning_data %>% 
    filter(estimate == 'mean') %>% 
    dplyr::select(-c(estimate, month_input, months)) %>%
    rowwise() %>% 
    mutate(ve_pred = if_else(prior_inf == 1, 
                             max(ve_pred - pe_nonsevere_level_hybrid - (weeks * pe_nonsevere_rate_hybrid/104), 0),
                             max(ve_pred - pe_nonsevere_level_vacc - (weeks * pe_nonsevere_rate_vacc/104), 0))) %>%
    rename(nonsevere_ve_pred = ve_pred)
  
  nonsevere_waning_prior_inf_only <- nonsevere_waning_prior_inf_only_data %>%
    filter(estimate == "mean") %>%
    rowwise() %>%
    mutate(ve_pred = max(ve_pred - pe_nonsevere_level_inf - (weeks * pe_nonsevere_rate_inf/104), 0)) %>%
    rename(nonsevere_prior_inf_only_ve_pred = ve_pred) %>%
    dplyr::select(-c(estimate, month_input, months, study))
  
  waning_data_clean <- merge(merge(setDT(severe_waning), setDT(nonsevere_waning), by = c("age_group", "prior_inf", "immunocompromised", "weeks"), all.x = TRUE),
                             merge(setDT(severe_waning_prior_inf_only), setDT(nonsevere_waning_prior_inf_only), by = c("weeks", "age_group", "immunocompromised"), all.x = TRUE),
                             by = c("weeks", "age_group", "immunocompromised"), all.x = TRUE)
  
  #Assigning protection to population
  entire_pop_with_pe <- merge(setDT(entire_pop), waning_data_clean, by.x = c("age_group", "prior_inf", "immunocompromised", "weeks_since_last_dose_inf"), by.y = c("age_group", "prior_inf", "immunocompromised", "weeks"), all.x = TRUE)
  
  #Adjust protection
  #Individuals who are unvaccinated and no prior infection history has no protection
  immune_naive_index <- which(entire_pop_with_pe$num_doses=="unvax" & entire_pop_with_pe$prior_inf == 0)
  entire_pop_with_pe[immune_naive_index, c("severe_ve_pred", "nonsevere_ve_pred")] <- 0
  
  #Individuals who are unvaccinated and have prior infection history have prior infection only waning immunity
  prior_inf_only_index <- which(entire_pop_with_pe$num_doses=="unvax" & entire_pop_with_pe$prior_inf == 1)
  entire_pop_with_pe[prior_inf_only_index, c("severe_ve_pred", "nonsevere_ve_pred")] <- entire_pop_with_pe[prior_inf_only_index, c("severe_prior_inf_only_ve_pred", "nonsevere_prior_inf_only_ve_pred")]
  
  #Individuals with infection <3 months from simulation start has perfect immunity
  perfect_immunity_index <- which(entire_pop_with_pe$weeks_since_last_inf<13 | entire_pop_with_pe$weeks_since_last_reinf<13)
  entire_pop_with_pe[perfect_immunity_index, c("severe_ve_pred", "nonsevere_ve_pred")] <- 1
  
  pe_by_age <- entire_pop_with_pe %>%
    group_by(age_group) %>%
    summarise(nonsevere_ve_pred = mean(nonsevere_ve_pred))
  
  #Setting average weekly incidence at start of simulation (averaged weekly estimates over Dec 2022)
  average_weekly_incidence <- average_severe_incidence %>%
    mutate(age_group = if_else(age_group == "0-17 years (Children)", "0-17 years", age_group)) %>%
    filter(weeks_since < 0,
           age_group != "All") %>%
    group_by(age_group) %>% summarise(avg_severe_inc = mean(inc)/100000) %>%
    mutate(avg_nonsevere_inc = c(avg_severe_inc[2], avg_severe_inc[2:5]) * c(200*1.25, 200*1.25, 79.6 *1.25, 22.6*1.25, 9.6*1.25) * nonsevere_mult,
           pe_avg = pe_by_age$nonsevere_ve_pred, 
           pe_minus_1 = 1 - pe_avg,
           pe_adj = (1 - pe_by_age$nonsevere_ve_pred[4])/pe_minus_1,
           inc_adj = avg_nonsevere_inc/avg_nonsevere_inc[4],
           beta = pe_adj * inc_adj,
           avg_severe_inc = if_else(age_group == "0-17 years", 0, avg_severe_inc))
  
  #Calculating total infections (dynamic term) at t = 0
  #NOTE: These are back-calculated from the severe/nonsevere incidences at t = 0
  inf_by_age <- merge(entire_pop_with_pe %>% group_by(age_group, immunocompromised) %>% summarise(total_pop = n()),
                      average_weekly_incidence, by = "age_group", all.x = TRUE) %>% 
    mutate(avg_severe_inc = if_else(immunocompromised %in% c(1, 2), avg_severe_inc * 2.8, avg_severe_inc),
           nonsevere_inf = total_pop * avg_nonsevere_inc,
           severe_inf = total_pop * avg_severe_inc,
           total_inf = nonsevere_inf + severe_inf) %>%
    group_by(age_group) %>% summarise(total_inf = ceiling(sum(total_inf)),
                                      total_pop = sum(total_pop))
  
  #Relative adjustment terms for contact matrix mixing
  contact_matrix_adj <- data.table(age_group = c("0-17 years","18-49 years", "50-64 years", "65-74 years", "75+ years"),
                                   contact_matrix_adj = (sum(inf_by_age$total_inf)/5000000)/c(sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X0.17.years),
                                                                                              sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X18.49.years),
                                                                                              sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X50.64.years),
                                                                                              sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X65.74.years),
                                                                                              sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X75..years)))
  
  #Multiplier adjustment terms 
  severe_infection_multipliers <- setDT(entire_pop_with_pe %>% 
                                          group_by(age_group, immunocompromised) %>% 
                                          summarise(mean_severe_ve = mean(severe_ve_pred),
                                                    mean_nonsevere_ve = mean(nonsevere_ve_pred)) %>%
                                          mutate(multiplier_adj = (1-mean_severe_ve)/(1-mean_nonsevere_ve)))
  severe_infection_multipliers$multiplier <- c(0, 1/c(200 *1.25, 79.6 *1.25, 22.6*1.25, 9.6*1.25))
  
  
  #Merge adjustment terms to waning data clean
  waning_data_clean <- merge(merge(merge(merge(waning_data_clean, severe_infection_multipliers, by = c("age_group", "immunocompromised"), all.x = TRUE),
                                   contact_matrix_adj, by = c("age_group"), all.x = TRUE),
                             setDT(average_weekly_incidence %>% dplyr::select(age_group, beta)), by = c("age_group"), all.x = TRUE),
                             setDT(hosp_death_age_stratified), by = c("age_group"), all.x = TRUE)
  
  
  #Assign vaccinations (for bivalent + XBB.1.5 to population)
  set.seed(88)
  vax_assigned <- realistic_vax_assignment(entire_pop_with_pe)
  
  rm(entire_pop_with_pe)
  
  return(list(waning_data_clean, vax_assigned, inf_by_age))
  
}


# 
# start.time <- Sys.time()
# #...Relevent codes...
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken    
