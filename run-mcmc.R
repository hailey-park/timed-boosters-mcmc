########################################################################################################################
#Title: MCMC Calibration
#Author: Hailey Park
#Date: April 16th, 2024
########################################################################################################################

rm(list=ls())
gc()

#Load libraries
library(dplyr)
library(ggplot2)
library(tibble)
library(lubridate)
library(here)
library(data.table)
library(tidyverse)

# Observed data
observed_data <- read.csv("data/clean-data/observed-inc-estimates.csv")[,-1] %>%
  filter(weeks_since > 0,
         age_group %in% c("18-49 years", "50-64 years", "65-74 years", "75+ years")) %>%
  dplyr::select(-c("week")) %>%
  group_by(age_group) %>%
  pivot_wider(names_from = age_group, values_from = inc)
  
#Set initial params
starting_params <- c(lambda = 1.13,
                     pe_nonsevere_level_vacc = 0,
                     pe_nonsevere_level_inf = 0,
                     pe_nonsevere_level_hybrid = 0,
                     pe_nonsevere_rate_vacc = 0,
                     pe_nonsevere_rate_inf = 0,
                     pe_nonsevere_rate_hybrid = 0,
                     nonsevere_mult = 2.5)

# Reading model scripts
source(here::here("model-setup-onetime.R"))
source(here::here("model-setup-eachtime.R"))
source(here::here("model-functions.R"))
source(here::here("mcmc-calibration-functions.R"))


#Set simulation size
sim_size <- 10000

#Run MCMC
chain = run_metropolis_MCMC(starting_params, sim_size)

burnIn = 1000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

#Set up folder structure to save simulation results
dir.create("simulation-results")

#Save results
write.csv(chain, "simulation-results/tester-mcmc-chain.csv")

# params["beta"] <- chain[sim_size, 2]
# params["gamma"] <- chain[sim_size, 1]
# times <- c(data$week)
# model.pred <- prediction(params,times)*params["N"]
# plot(incidence~week,data=data,type='b',col='red') 
# lines(times,model.pred,type='l')



