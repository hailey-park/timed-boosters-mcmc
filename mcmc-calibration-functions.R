########################################################################################################################
#Title: MCMC Calibration
#Author: Hailey Park
#Date: April 16th, 2024
########################################################################################################################


prediction <- function (params) {

  #Set varied params (lambda, pe vars, nonsevere mult, etc.)
  lambda <- params[1]
  pe_nonsevere_level_vacc <- params[2]
  pe_nonsevere_level_inf <- params[3]
  pe_nonsevere_level_hybrid <- params[4]
  pe_nonsevere_rate_vacc <- params[5]
  pe_nonsevere_rate_inf <- params[6]
  pe_nonsevere_rate_hybrid <- params[7]
  nonsevere_mult <- params[8]

  #Adjust waning curves and other adjustments to model
  source(here::here("model-setup-eachtime.R"))
  
  #Run model
  set.seed(88)
  vax_assigned <- realistic_vax_assignment(entire_pop)
  
  set.seed(88)
  results <- boosterSimulation(vax_assigned %>% mutate(lambda = lambda))
  
  reformatted_results <- results
  
  sim_results <- melt(colSums(sim_data_raw %>% select(week1:week70)))
  sim_results$weeks <- rownames(sim_results)
  sim_results <- sim_results %>%
    mutate(weeks = as.numeric(str_sub(weeks, 5)) + 12,
           simulated_inc = value/10000000 * 100000) %>%
    add_row(weeks = 0,
            simulated_inc = NA)
  
  
  return(results)
}


# Define likelihood function with poisson distribution (count data)
poisson.loglik <- function (params, data) {

  pred <- prediction(params)
  pred_I2 <- matrix(nrow=dim(pred)[1], ncol=dim(pred)[2])
  pred_I2[1,] <- rep(0,dim(pred)[2])
  for (j in 2:dim(pred)[1]) {
    pred_I2[j,] <- as.numeric(pred[j,] - pred[j-1,])
  }
  pred_I <- matrix(nrow=(dim(pred)[1]/12), ncol=dim(pred)[2])
  for (j in 1:(dim(pred)[1]/12)) {
    pred_I[j,] <- colSums(pred_I2[((j-1)*12+1):(j*12),])
  }
  pred_I[1,] <- pred_I[1,]/(11/12)
  # For MCMC
  # Use subset of age groups (use >1 case per month)
  # adjust for population (0.95 correction)
  age_18_49 <- sum(dpois(x=(round(tail(data$incidence1*12, 10))+1),lambda=((pred_I[,1])+1), log=TRUE))
  age_50_64 <- sum(dpois(x=(round(tail(data$incidence2*12, 10))+1),lambda=((pred_I[,2])+1), log=TRUE))
  age_65_74 <- sum(dpois(x=(round(tail(data$incidence3*12, 10))+1),lambda=((pred_I[,3])+1), log=TRUE))
  age_75_plus <- sum(dpois(x=(round(tail(data$incidence4*12, 10))+1),lambda=((pred_I[,4])+1), log=TRUE))

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
  
  lambda_prior = dnorm(lambda, mean=0.9, sd=0.12, log = T)
  # pe_severe_level_vacc_prior = dnorm(pe_severe_level_vacc, mean=0, sd=0.1, log = F)
  # pe_severe_level_inf_prior = dnorm(pe_severe_level_inf, mean=0, sd=0.3, log = F)
  # pe_severe_level_hybrid_prior = dnorm(pe_severe_level_hybrid, mean=0.025, sd=0.03, log = F)
  # pe_severe_rate_vacc_prior = dnorm(pe_severe_rate_vacc, mean=0, sd=0.1, log = F)
  # pe_severe_rate_inf_prior = dnorm(pe_severe_rate_inf, mean=0.025, sd=0.03, log = F)
  # pe_severe_rate_hybrid_prior = dnorm(pe_severe_rate_hybrid, mean=0, sd=0.2, log = F)
  pe_nonsevere_level_vacc_prior = dnorm(pe_nonsevere_level_vacc, mean=0, sd=0.1, log = F)
  pe_nonsevere_level_inf_prior = dnorm(pe_nonsevere_level_inf, mean=0, sd=0.1, log = F)
  pe_nonsevere_level_hybrid_prior = dnorm(pe_nonsevere_level_hybrid, mean=0, sd=0.1, log = F)
  pe_nonsevere_rate_vacc_prior = dnorm(pe_nonsevere_rate_vacc, mean=0, sd=0.1, log = F)
  pe_nonsevere_rate_inf_prior = dnorm(pe_nonsevere_rate_inf, mean=0, sd=0.04, log = F)
  pe_nonsevere_rate_hybrid_prior = dnorm(pe_nonsevere_rate_hybrid, mean=0, sd=0.06, log = F)
  nonsevere_mult_prior = dnorm(nonsevere_mult, mean=0.9, sd=0.15, log = T)
  
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
  return(rnorm(8,mean = param, sd= c(0.1, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.1)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}
