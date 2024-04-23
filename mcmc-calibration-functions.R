########################################################################################################################
#Title: MCMC Calibration
#Author: Hailey Park
#Date: April 16th, 2024
########################################################################################################################

prediction <- function (params) {

  #print(params)
  #Set varied params (lambda, pe vars, nonsevere mult, etc.)
  # lambda <- params[1]
  # pe_nonsevere_level_vacc <- params[2]
  # pe_nonsevere_level_inf <- params[3]
  # pe_nonsevere_level_hybrid <- params[4]
  # pe_nonsevere_rate_vacc <- params[5]
  # pe_nonsevere_rate_inf <- params[6]
  # pe_nonsevere_rate_hybrid <- params[7]
  # nonsevere_mult <- params[8]

  #Adjust waning curves and other adjustments to model
  model_set_up <- model_set_up_fn(params)
  waning_data_clean <- model_set_up[[1]]
  vax_assigned <- model_set_up[[2]]
  inf_by_age <- model_set_up[[3]]

  #Run model
  #set.seed(88)
  results <- boosterSimulation(vax_assigned %>% mutate(lambda = params[1]), waning_data_clean, inf_by_age) 
  
  #Reformat results into weekly incidence estimates by age group (separate columns)
  reformatted_results <- melt(setDT(results %>% select(total_pop, week1:week52)), id = c("age_group", "total_pop")) %>%
    mutate(weeks = as.numeric(str_sub(variable, 5)),
           simulated_inc = value/total_pop * 100000) %>%
    dplyr::select(-c(variable, value, total_pop)) %>%
    group_by(age_group) %>%
    pivot_wider(names_from = age_group, values_from = simulated_inc)

    return(reformatted_results)
}


# Define likelihood function with poisson distribution (count data)
poisson.loglik <- function (params, data) {

  sim_pred <- prediction(params)

    # For MCMC
  # Use subset of age groups (use >1 case per month)
  age_18_49 <- sum(dpois(x=(round(observed_data$`18-49 years`)+1),lambda=(round(sim_pred$`18-49 years`)+1), log=TRUE))
  age_50_64 <- sum(dpois(x=(round(observed_data$`50-64 years`)+1),lambda=(round(sim_pred$`50-64 years`)+1), log=TRUE))
  age_65_74 <- sum(dpois(x=(round(observed_data$`65-74 years`)+1),lambda=(round(sim_pred$`65-74 years`)+1), log=TRUE))
  age_75_plus <- sum(dpois(x=(round(observed_data$`75+ years`)+1),lambda=(round(sim_pred$`75+ years`)+1), log=TRUE))

  total <- sum(age_18_49, age_50_64, age_65_74, age_75_plus)
  return(total)
}


#The likelihood is the probability (density) with which we would expect the 
# observed data to occur conditional on the parameters of the model that we look at.
likelihood <- function (par_initial) {
  par <- par_initial
  poisson.loglik(par,observed_data)
}


# Prior distribution
prior <- function(par_initial){
  lambda = par_initial[1]
  # pe_severe_level_vacc = par_initial[2] #Don't use these for now
  # pe_severe_level_inf = par_initial[3]
  # pe_severe_level_hybrid = par_initial[4]
  # pe_severe_rate_vacc = par_initial[5]
  # pe_severe_rate_inf = par_initial[6]
  # pe_severe_rate_hybrid = par_initial[7]
  pe_nonsevere_level_vacc = par_initial[2]
  pe_nonsevere_level_inf = par_initial[3]
  pe_nonsevere_level_hybrid = par_initial[4]
  pe_nonsevere_rate_vacc = par_initial[5]
  pe_nonsevere_rate_inf = par_initial[6]
  pe_nonsevere_rate_hybrid = par_initial[7]
  nonsevere_mult = par_initial[8]

  lambda_prior = dlnorm(lambda, mean=0.9, sd=0.12, log = T)
  # pe_severe_level_vacc_prior = dnorm(pe_severe_level_vacc, mean=0, sd=0.1, log = F)
  # pe_severe_level_inf_prior = dnorm(pe_severe_level_inf, mean=0, sd=0.3, log = F)
  # pe_severe_level_hybrid_prior = dnorm(pe_severe_level_hybrid, mean=0.025, sd=0.03, log = F)
  # pe_severe_rate_vacc_prior = dnorm(pe_severe_rate_vacc, mean=0, sd=0.1, log = F)
  # pe_severe_rate_inf_prior = dnorm(pe_severe_rate_inf, mean=0.025, sd=0.03, log = F)
  # pe_severe_rate_hybrid_prior = dnorm(pe_severe_rate_hybrid, mean=0, sd=0.2, log = F)
  pe_nonsevere_level_vacc_prior = dnorm(pe_nonsevere_level_vacc, mean=0, sd=0.1, log = T)
  pe_nonsevere_level_inf_prior = dnorm(pe_nonsevere_level_inf, mean=0, sd=0.1, log = T)
  pe_nonsevere_level_hybrid_prior = dnorm(pe_nonsevere_level_hybrid, mean=0, sd=0.1, log = T)
  pe_nonsevere_rate_vacc_prior = dnorm(pe_nonsevere_rate_vacc, mean=0, sd=0.1, log = T)
  pe_nonsevere_rate_inf_prior = dnorm(pe_nonsevere_rate_inf, mean=0, sd=0.1, log = T)
  pe_nonsevere_rate_hybrid_prior = dnorm(pe_nonsevere_rate_hybrid, mean=0, sd=0.1, log = T)
  nonsevere_mult_prior = dlnorm(nonsevere_mult, mean=0.9, sd=0.15, log = T)
  
  return(lambda_prior+
           # pe_severe_level_vacc_prior+pe_severe_level_inf_prior+pe_severe_level_hybrid_prior+
           # pe_severe_rate_vacc_prior+pe_severe_rate_inf_prior+pe_severe_rate_hybrid_prior+
           pe_nonsevere_level_vacc_prior+pe_nonsevere_level_inf_prior+pe_nonsevere_level_hybrid_prior+
           pe_nonsevere_rate_vacc_prior+pe_nonsevere_rate_inf_prior+pe_nonsevere_rate_hybrid_prior+
           nonsevere_mult_prior)
}


posterior <- function(param){
  return (likelihood(param) + prior(param))
}


#Choosing a new parameter value close to the old value based on some 
#probability density that is called the proposal function
proposalfunction <- function(param){
  return(rnorm(8,mean = param, sd= c(0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.2)))
}


run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,8))
  chain[1,] = startvalue
  print(startvalue)
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    print(proposal)
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    
    print(chain[i+1,])
  }
  return(chain)
}
