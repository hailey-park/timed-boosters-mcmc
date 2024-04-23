########################################################################################################################
#Title: Model Functions
#Author: Hailey Park
#Date: April 19th, 2024
########################################################################################################################


#Function for outcome occurrence based on risk (Risk = Lambda* (1 - PE))
outcome_occurrence <- function(age, inf, time, immuno, doses, lambda, perfect_immunity_counter, death_marker, weekly_infection_by_age, contact_matrix_adj_factors, waning_data_clean, inf_by_age) {
  
  severe_pe <- rep(1, length(age))
  nonsevere_pe <- rep(1, length(age))
  beta <- rep(0, length(age))
  severe_multiplier_with_adj <- rep(1, length(age))
  
  #Creating a df of individuals eligible for infection to merge with waning_data_clean to get protection at specific time point
  index_individuals_eligible <- which(perfect_immunity_counter == 0 & death_marker == 0)
  df_individuals_eligible <- data.table(index_individual = index_individuals_eligible,
                                        age_group = age[index_individuals_eligible],
                                        prior_inf = inf[index_individuals_eligible],
                                        immunocompromised = immuno[index_individuals_eligible],
                                        num_doses = doses[index_individuals_eligible],
                                        weeks = time[index_individuals_eligible],
                                        key = c("age_group", "prior_inf", "immunocompromised", "weeks"))
  
  # df_protection <- (df_individuals_eligible[waning_data_clean,
  #                                           on=c("age_group", "prior_inf", "immunocompromised", "weeks"),
  #                                           nomatch = NULL]) %>% arrange(index_individual)
  
  df_protection <- as.data.frame(merge(df_individuals_eligible, waning_data_clean, by=c("age_group", "prior_inf", "immunocompromised", "weeks"), all.x = TRUE)) %>%
    arrange(index_individual)
  

  # print(df_protection %>% filter(index_individual == 500000))
  
  #For individuals who are immune naive, they have 0 protection
  df_protection$severe_ve_pred[df_protection$num_doses == "unvax" & df_protection$prior_inf == 0] <- 0
  df_protection$nonsevere_ve_pred[df_protection$num_doses == "unvax" & df_protection$prior_inf == 0] <- 0
  
  #Individuals who are unvaccinated and have prior infection history have prior infection only waning immuntiy
  df_protection$severe_ve_pred[df_protection$num_doses=="unvax" & df_protection$prior_inf == 1] <- df_protection$severe_prior_inf_only_ve_pred[df_protection$num_doses=="unvax" & df_protection$prior_inf == 1]
  df_protection$nonsevere_ve_pred[df_protection$num_doses=="unvax" & df_protection$prior_inf == 1] <- df_protection$nonsevere_prior_inf_only_ve_pred[df_protection$num_doses=="unvax" & df_protection$prior_inf == 1]

  #Calculate protection (severe + nonsevere), beta, and multipliers for individuals eligible for infection
  severe_pe[index_individuals_eligible] <- df_protection$severe_ve_pred
  nonsevere_pe[index_individuals_eligible] <- df_protection$nonsevere_ve_pred
  beta[index_individuals_eligible] <- df_protection$beta
  severe_multiplier_with_adj[index_individuals_eligible] <- df_protection$multiplier/df_protection$multiplier_adj
  
  #Calculate the age contact matrix terms for all individuals
  contact_matrix_term_by_age <- data.table(age_group = c("0-17 years","18-49 years", "50-64 years", "65-74 years", "75+ years"),
                                           contact_matrix_term = c(sum(weekly_infection_by_age/inf_by_age$total_pop * contact_matrix$X0.17.years),
                                                                    sum(weekly_infection_by_age/inf_by_age$total_pop * contact_matrix$X18.49.years),
                                                                    sum(weekly_infection_by_age/inf_by_age$total_pop * contact_matrix$X50.64.years),
                                                                    sum(weekly_infection_by_age/inf_by_age$total_pop * contact_matrix$X65.74.years),
                                                                    sum(weekly_infection_by_age/inf_by_age$total_pop * contact_matrix$X75..years)))
  
  weekly_contact_matrix_term <- (data.table(age_group = age)[contact_matrix_term_by_age, 
                                                        on=c("age_group"), 
                                                        nomatch = NULL])$contact_matrix_term

  nonsevere_risk <- lambda * (1 - nonsevere_pe) * beta * contact_matrix_adj_factors * (weekly_contact_matrix_term)
  
  severe_risk <- lambda * (1 - severe_pe) * beta * contact_matrix_adj_factors * (weekly_contact_matrix_term) * severe_multiplier_with_adj
  
  #if nonsevere or severe risk > 1, set to 1
  nonsevere_risk[nonsevere_risk > 1] <- 1
  severe_risk[severe_risk > 1] <- 1
  
  #Simulate outcomes
  severe_outcomes <- rbinom(length(severe_risk), 1, severe_risk)
  nonsevere_outcomes <- rbinom(length(nonsevere_risk), 1, nonsevere_risk)
  
  return(list(severe_outcomes, nonsevere_outcomes))
}

########################################################################################################################
#Vaccination Assignment under current coverages
#Assign who is getting vaccinated using age-specific realistic coverage rates. This is important for the 1-booster base case
#Timing of vaccination is assigned inside the vaccination intervention functions.

realistic_vax_assignment <- function(df){
  #Assign who is getting vaccinated using age-specific coverage rates
  vax_assignment <- df %>% mutate(biv_vax = 0,
                                  xbb_vax = 0)
  
  age_0_17_index <- which(vax_assignment$age_group == "0-17 years" & vax_assignment$num_doses != "unvax")
  age_18_49_index <- which(vax_assignment$age_group == "18-49 years" & vax_assignment$num_doses != "unvax")
  age_50_64_index <- which(vax_assignment$age_group == "50-64 years" & vax_assignment$num_doses != "unvax")
  age_65_plus_index <- which((vax_assignment$age_group == "65-74 years" | vax_assignment$age_group == "75+ years") & vax_assignment$num_doses != "unvax")
  immunocompromised_index <- which(vax_assignment$immunocompromised %in% c(1, 2)  & vax_assignment$num_doses != "unvax")
  
  #Bivalent coverage maximum uptake by age groups
  vax_assignment$biv_vax[age_0_17_index] <- rbinom(length(age_0_17_index), 1, 0.1757 * length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$num_doses != "unvax" & vax_assignment$immunocompromised == 0)) / length(age_0_17_index))
  vax_assignment$biv_vax[age_18_49_index] <- rbinom(length(age_18_49_index), 1, 0.2103 * length(which(vax_assignment$age_group == "18-49 years" & vax_assignment$num_doses != "unvax" & vax_assignment$immunocompromised == 0)) / length(age_18_49_index))
  vax_assignment$biv_vax[age_50_64_index] <- rbinom(length(age_50_64_index), 1, 0.334 * length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$num_doses != "unvax" & vax_assignment$immunocompromised == 0)) / length(age_50_64_index))
  vax_assignment$biv_vax[age_65_plus_index] <- rbinom(length(age_65_plus_index), 1, 0.5413 )#* length(which(vax_assignment$age_group %in% c("65-74 years", "75+ years") & vax_assignment$num_doses != "unvax" & vax_assignment$immunocompromised == 0)) / length(age_65_plus_index))
  vax_assignment$biv_vax[immunocompromised_index] <- rbinom(length(immunocompromised_index), 1, 0.6) #fake assumption; revisit
  
  # print(paste0("Percentage of 0-17 years receiving bivalent vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$num_doses != "unvax" & vax_assignment$biv_vax == 1))/length(age_0_17_index)))
  # print(paste0("Percentage of 18-49 years receiving bivalent vaccinated: ", length(which(vax_assignment$age_group == "18-49 years" & vax_assignment$num_doses != "unvax" & vax_assignment$biv_vax == 1))/length(age_18_49_index)))
  # print(paste0("Percentage of 50-64 years receiving bivalent vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$num_doses != "unvax" & vax_assignment$biv_vax == 1))/length(age_50_64_index)))
  # print(paste0("Percentage of 65+ years receiving bivalent vaccinated: ",  length(which(vax_assignment$age_group %in% c("65-74 years", "75+ years") & vax_assignment$num_doses != "unvax" & vax_assignment$biv_vax == 1))/length(age_65_plus_index)))
  # print(paste0("Percentage of immunocompromised receiving bivalent vaccinated: ",  length(which(vax_assignment$immunocompromised %in% c(1, 2) & vax_assignment$num_doses != "unvax" & vax_assignment$biv_vax == 1))/length(immunocompromised_index)))

  #XBB.15 coverage maximum uptake by age groups (only administering XBB booster to those who already received bivalent booster)
  vax_assignment$xbb_vax[intersect(age_0_17_index, which(vax_assignment$biv_vax == 1))] <- rbinom(length(intersect(age_0_17_index, which(vax_assignment$biv_vax == 1))), 1, 0.0957 * length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$num_doses != "unvax" & vax_assignment$biv_vax == 1 & vax_assignment$immunocompromised == 0)) / length(intersect(age_0_17_index, which(vax_assignment$biv_vax == 1))) * (length(age_0_17_index)/length(intersect(age_0_17_index, which(vax_assignment$biv_vax == 1)))))
  vax_assignment$xbb_vax[intersect(age_18_49_index, which(vax_assignment$biv_vax == 1))] <- rbinom(length(intersect(age_18_49_index, which(vax_assignment$biv_vax == 1))), 1, 0.0974 * length(which(vax_assignment$age_group == "18-49 years" & vax_assignment$num_doses != "unvax" & vax_assignment$biv_vax == 1 & vax_assignment$immunocompromised == 0)) / length(intersect(age_18_49_index, which(vax_assignment$biv_vax == 1))) * (length(age_18_49_index)/length(intersect(age_18_49_index, which(vax_assignment$biv_vax == 1)))))
  vax_assignment$xbb_vax[intersect(age_50_64_index, which(vax_assignment$biv_vax == 1))] <- rbinom(length(intersect(age_50_64_index, which(vax_assignment$biv_vax == 1))), 1, 0.1807 * length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$num_doses != "unvax" & vax_assignment$biv_vax == 1 & vax_assignment$immunocompromised == 0)) / length(intersect(age_50_64_index, which(vax_assignment$biv_vax == 1))) * (length(age_50_64_index)/length(intersect(age_50_64_index, which(vax_assignment$biv_vax == 1)))))
  vax_assignment$xbb_vax[intersect(age_65_plus_index, which(vax_assignment$biv_vax == 1))] <- rbinom(length(intersect(age_65_plus_index, which(vax_assignment$biv_vax == 1))), 1, 0.347  * (length(age_65_plus_index)/length(intersect(age_65_plus_index, which(vax_assignment$biv_vax == 1))))) #* length(which(vax_assignment$age_group %in% c("65-74 years", "75+ years") & vax_assignment$num_doses != "unvax" & vax_assignment$biv_vax == 1 & vax_assignment$immunocompromised == 0)) / length(intersect(age_65_plus_index, which(vax_assignment$biv_vax == 1)))
  vax_assignment$xbb_vax[intersect(immunocompromised_index, which(vax_assignment$biv_vax == 1))] <- rbinom(length(intersect(immunocompromised_index, which(vax_assignment$biv_vax == 1))), 1, 0.4 * (length(immunocompromised_index)/length(intersect(immunocompromised_index, which(vax_assignment$biv_vax == 1))))) #fake assumption; revisit
  
  # print(paste0("Percentage of 0-17 years receiving XBB vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$num_doses != "unvax" & vax_assignment$xbb_vax == 1))/length(age_0_17_index)))
  # print(paste0("Percentage of 18-49 years receiving XBB vaccinated: ", length(which(vax_assignment$age_group == "18-49 years" & vax_assignment$num_doses != "unvax" & vax_assignment$xbb_vax == 1))/length(age_18_49_index)))
  # print(paste0("Percentage of 50-64 years receiving XBB vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$num_doses != "unvax" & vax_assignment$xbb_vax == 1))/length(age_50_64_index)))
  # print(paste0("Percentage of 65+ years receiving XBB vaccinated: ",  length(which(vax_assignment$age_group %in% c("65-74 years", "75+ years") & vax_assignment$num_doses != "unvax" & vax_assignment$xbb_vax == 1))/length(age_65_plus_index)))
  # print(paste0("Percentage of immunocompromised receiving XBB vaccinated: ",  length(which(vax_assignment$immunocompromised %in% c(1, 2) & vax_assignment$num_doses != "unvax" & vax_assignment$xbb_vax == 1))/length(immunocompromised_index)))
  # 
  return(vax_assignment %>% arrange(individual))
}

########################################################################################################################
boosterSimulation <- function(df, waning_data_clean, inf_by_age){
  
  
  #Store severe and nonsevere outcome counts in df
  grouped_outcome_counts <- df  %>% group_by(age_group, immunocompromised) %>% summarise(total_pop = n())
  grouped_outcome_counts[sprintf("week%s",(1:52))] <- NA
  grouped_outcome_counts[sprintf("nonsevere_week%s",(1:52))] <- NA
  
  #Input data (entire pop)
  input <- df %>% arrange(individual)
  input[sprintf("week%s",(1))] <- NA
  input[sprintf("nonsevere_week%s",(1))] <- NA
  
  #merge with waning data clean to get other factors
  df <-  as.data.frame(merge(setDT(df %>% dplyr::select(individual, age_group)), 
                             waning_data_clean %>% group_by(age_group) %>% summarise(contact_matrix_adj = mean(contact_matrix_adj),
                                                                                     perc_death = mean(perc_death)),
                             by=c("age_group"), 
                             all.x = TRUE)) %>%
    arrange(individual)
  
  #Population's info (age_group, num_doses, prior_inf, etc.) at each timestep
  age <- as.character(input$age_group)
  doses <- as.character(input$num_doses)
  inf <- input$prior_inf
  immuno <- input$immunocompromised
  doses <- input$num_doses
  time_since_last <- input$weeks_since_last_dose_inf
  time_since_last_dose <- as.numeric(as.character(input$weeks_since_last_dose))
  lambda <- input$lambda
  # contact_matrix_adj_factors <-(data.table(age_group = age)[contact_matrix_adj,
  #                                                           on=c("age_group"), 
  #                                                           nomatch = NULL])$contact_matrix_adj
  # 
  # prob_death <- (data.table(age_group = age)[hosp_death_age_stratified,
  #                                            on=c("age_group"), 
  #                                            nomatch = NULL])$perc_death
  
  contact_matrix_adj_factors <-df$contact_matrix_adj
  
  prob_death <- df$perc_death
  
  perfect_immunity_counter <- rep(0,nrow(input)) #If non-death infection occurs, counting down perfect immunity weeks
  index_recent_infection <- which(inf == 1 & time_since_last < 13 & time_since_last < time_since_last_dose) #Individuals infected in 3 months preceding start of sim have perfect immunity at start
  perfect_immunity_counter[index_recent_infection] <- 14 - time_since_last[index_recent_infection] #REVISIT
  death_marker <- rep(0,nrow(input)) #If death occurs
  hosp_count <- rep(0, nrow(input))
  death_count <- rep(0, nrow(input))
  weeks <- c(1:104)
  
  #get indexes by age group and by bivalent/xbb vaccine coverages
  to_vaccinate_biv_index <- which(input$biv_vax == 1)
  to_vaccinate_xbb_index <- which(input$xbb_vax == 1)
  age_0_17_index <- which(input$age_group == "0-17 years")
  age_18_49_index <- which(input$age_group == "18-49 years")
  age_50_64_index <- which(input$age_group == "50-64 years")
  age_65_plus_index <- which((input$age_group == "65-74 years" | input$age_group == "75+ years"))
  
  #distribute bivalent and xbb booster uptake over 37 weeks (bivalent) and 16 weeks (XBB), according to age distributions
  vaccine_wave_biv <- rep(0, nrow(input))
  vaccine_wave_xbb <- rep(0, nrow(input))
  vaccine_wave_biv[intersect(to_vaccinate_biv_index, age_0_17_index)] <- sample(c(1:37), length(intersect(to_vaccinate_biv_index, age_0_17_index)), prob = (bivalent_coverage %>% filter(age_categories == "0-17") %>% arrange(week))$prop_of_perc_up_to_date, replace = TRUE)
  vaccine_wave_biv[intersect(to_vaccinate_biv_index, age_18_49_index)] <- sample(c(1:37), length(intersect(to_vaccinate_biv_index, age_18_49_index)), prob = (bivalent_coverage %>% filter(age_categories == "18-49") %>% arrange(week))$prop_of_perc_up_to_date, replace = TRUE)
  vaccine_wave_biv[intersect(to_vaccinate_biv_index, age_50_64_index)] <- sample(c(1:37), length(intersect(to_vaccinate_biv_index, age_50_64_index)), prob = (bivalent_coverage %>% filter(age_categories == "50-64") %>% arrange(week))$prop_of_perc_up_to_date, replace = TRUE)
  vaccine_wave_biv[intersect(to_vaccinate_biv_index, age_65_plus_index)] <- sample(c(1:37), length(intersect(to_vaccinate_biv_index, age_65_plus_index)), prob = (bivalent_coverage %>% filter(age_categories == "65+") %>% arrange(week))$prop_of_perc_up_to_date, replace = TRUE)
  
  vaccine_wave_xbb[intersect(to_vaccinate_xbb_index, age_0_17_index)] <- sample(c(38:53), length(intersect(to_vaccinate_xbb_index, age_0_17_index)), prob = (xbb_coverage %>% filter(age_categories == "0-17") %>% arrange(week))$prop_of_perc_up_to_date, replace = TRUE)
  vaccine_wave_xbb[intersect(to_vaccinate_xbb_index, age_18_49_index)] <- sample(c(38:53), length(intersect(to_vaccinate_xbb_index, age_18_49_index)), prob = (xbb_coverage %>% filter(age_categories == "18-49") %>% arrange(week))$prop_of_perc_up_to_date, replace = TRUE)
  vaccine_wave_xbb[intersect(to_vaccinate_xbb_index, age_50_64_index)] <- sample(c(38:53), length(intersect(to_vaccinate_xbb_index, age_50_64_index)), prob = (xbb_coverage %>% filter(age_categories == "50-64") %>% arrange(week))$prop_of_perc_up_to_date, replace = TRUE)
  vaccine_wave_xbb[intersect(to_vaccinate_xbb_index, age_65_plus_index)] <- sample(c(38:53), length(intersect(to_vaccinate_xbb_index, age_65_plus_index)), prob = (xbb_coverage %>% filter(age_categories == "65+") %>% arrange(week))$prop_of_perc_up_to_date, replace = TRUE)
  
  weekly_infection_by_age <- inf_by_age$total_inf

  #Iterate through each time step
  for (i in (1:52)) {
    # print(paste0("Week: ", i))
    
    #Staggering bivalent vaccination over 37-week window
    if(i %in% c(1:37)){
      vaccine_wave_index <- which(vaccine_wave_biv == i)
      time_since_last[vaccine_wave_index] <- 1
    }

    #Staggering XBB vaccination over 16-week window
    if(i %in% c(38:53)){
      vaccine_wave_index <- which(vaccine_wave_xbb == i)
      time_since_last[vaccine_wave_index] <- 1
    }

    time_since_last[time_since_last >= 104] <- 104     #Assuming that >24 month waning is same as 24 month waning pe
    time_since_last[time_since_last == 0] <- 1     
    
    week <- weeks[time_since_last]
    
    #Do outcomes occur?
    outcomes <- outcome_occurrence(age, inf, week, immuno, doses, lambda, perfect_immunity_counter, death_marker, weekly_infection_by_age, contact_matrix_adj_factors, waning_data_clean, inf_by_age)
    severe_outcomes <- outcomes[[1]]
    nonsevere_outcomes <- outcomes[[2]]
    
    # print(paste0("Total weekly severe infections: ", sum(severe_outcomes)))
    # print(paste0("Total weekly nonsevere infections: ", sum(nonsevere_outcomes)))
    
    #If no outcome occurs, increase time since last
    index_no_outcome <- which(severe_outcomes == 0 & nonsevere_outcomes == 0)
    time_since_last[index_no_outcome] <- time_since_last[index_no_outcome] + 1
    
    #Decrease 1 from perfect immunity counter (if applicable)
    perfect_immunity_counter[perfect_immunity_counter > 0] <- perfect_immunity_counter[perfect_immunity_counter > 0] - 1
    
    #If outcome occurs, 
    #change their prior infection status to 1, time since last to 1, perfect immunity counter to 13
    index_outcome <- which(severe_outcomes == 1 | nonsevere_outcomes == 1)
    inf[index_outcome] <- 1
    time_since_last[index_outcome] <- 1
    perfect_immunity_counter[index_outcome] <- 13
    
    #Then check if severe outcome is hosp vs. death
    index_severe_outcome <- which(severe_outcomes == 1)
    death_ind <- rbinom(length(index_severe_outcome), 1, prob_death[index_severe_outcome])
    
    #If death, cut simulation for individual (death marker)
    death_ind_index <- index_severe_outcome[which(death_ind == 1)]
    death_marker[death_ind_index] <- 1
    death_count[death_ind_index] <- death_count[death_ind_index] + 1
    
    hosp_ind_index <- index_severe_outcome[which(death_ind == 0)]
    hosp_count[hosp_ind_index] <- hosp_count[hosp_ind_index] + 1
    
    #If both severe outcome and nonsevere outcome occur in same individual, remove nonsevere outcome
    index_both_outcome <- which(severe_outcomes == 1 & nonsevere_outcomes == 1)
    nonsevere_outcomes[index_both_outcome] <- 0
    
    #Re-update weekly_infection_by_age counter with new infection counts
    weekly_infection_by_age[1] <- sum(severe_outcomes[which(age == "0-17 years")]) + sum(nonsevere_outcomes[which(age == "0-17 years")])
    weekly_infection_by_age[2] <- sum(severe_outcomes[which(age == "18-49 years")]) + sum(nonsevere_outcomes[which(age == "18-49 years")])
    weekly_infection_by_age[3] <- sum(severe_outcomes[which(age == "50-64 years")]) + sum(nonsevere_outcomes[which(age == "50-64 years")])
    weekly_infection_by_age[4] <- sum(severe_outcomes[which(age == "65-74 years")]) + sum(nonsevere_outcomes[which(age == "65-74 years")])
    weekly_infection_by_age[5] <- sum(severe_outcomes[which(age == "75+ years")]) + sum(nonsevere_outcomes[which(age == "75+ years")])
    
    #Add data to dataframe
    input$week1 <- severe_outcomes
    input$nonsevere_week1 <- nonsevere_outcomes
 
    grouped_outcomes <- input %>% group_by(age_group, immunocompromised) %>% summarise(total_severe = sum(week1),
                                                                                       total_nonsevere = sum(nonsevere_week1))
                                                                                       
    grouped_outcome_counts[, i + 3] <- grouped_outcomes$total_severe
    grouped_outcome_counts[, i + 55] <- grouped_outcomes$total_nonsevere

  }
  
  return(grouped_outcome_counts)
}

########################################################################################
# 
# # 
# set.seed(88)
# inspection <- clean_df %>%
#   lapply(realistic_vax_assignment)
# 
# 
# 
# set.seed(188)
# results <- boosterSimulation(inspection[[1]] %>% mutate(lambda = 1.13))
# write.csv(results, paste0("simulation-results/2.5x-updated sero-10% waning/lambda-1.03.csv"))
# 
# 
# results$week1 / results$total_pop
# 
# check <- read.csv("simulation-results/2.5x-updated sero/lambda-1.13-no vax.csv")[,-1]
# check$week1 / check$total_pop
