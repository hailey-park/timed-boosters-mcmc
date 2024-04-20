###################################################################################################
#Title: Model set-up (onetime script)
#Author: Hailey Park
#Date: April 17, 2024
###################################################################################################

#Load data 
pop_by_age <- read.csv("data/clean-data/pop_by_age.csv")[,-1]
cases_by_week <- read.csv("data/clean-data/cases_by_week.csv")[,-1]
four_doses_by_week <- read.csv("data/clean-data/four_doses_by_week.csv")[,-1]
three_doses_by_week <- read.csv("data/clean-data/three_doses_by_week.csv")[,-1]
contact_matrix <- read.csv("data/clean-data/contact matrix.csv")
hosp_death_age_stratified <- data.table(read.csv("data/clean-data/hosp_death_age_stratified_counts_adj.csv")[,-1]) 
vaccine_coverage_by_age <- read.csv("data/clean-data/vaccine-coverage-by-age.csv")[,-1]

#Clean data
four_doses_by_week$week <- as.character(four_doses_by_week$week)
three_doses_by_week$week <- as.character(three_doses_by_week$week)
contact_matrix <- contact_matrix %>% mutate(across(X0.17.years:X75..years) * 7)

#Vaccine coverage by age. calculating prop administered each week for simulations
bivalent_coverage <- vaccine_coverage_by_age %>% filter(administered_date <= as.Date("2023-09-13")) %>%
  mutate(week = 1 + as.numeric(as.Date(administered_date) - as.Date("2022-12-31")) %/% 7) %>%
  group_by(week, age_categories) %>%
  summarise(perc_up_to_date = max(perc_up_to_date)) %>%
  arrange(week) %>%
  group_by(age_categories) %>%
  mutate(diff = c((perc_up_to_date-lag(perc_up_to_date))[2], (perc_up_to_date-lag(perc_up_to_date))[-1]),
         total = sum(diff),
         prop_of_perc_up_to_date = diff/total * 100)

xbb_coverage <- vaccine_coverage_by_age %>% filter(administered_date > as.Date("2023-09-13")) %>%
  mutate(week = 38 + as.numeric(as.Date(administered_date) - as.Date("2023-09-14")) %/% 7) %>%
  group_by(week, age_categories) %>%
  summarise(perc_up_to_date = max(perc_up_to_date)) %>%
  arrange(week) %>%
  group_by(age_categories) %>%
  mutate(diff = c(perc_up_to_date[1], (perc_up_to_date-lag(perc_up_to_date))[-1]),
         total = sum(diff),
         prop_of_perc_up_to_date = diff/total * 100)

############################################################################
#If population df already exists, call it, otherwise, create it

if(file.exists("data/clean-data/entire_pop_5mil.csv")){
  entire_pop <- read.csv("data/clean-data/entire_pop_5mil.csv")[,-1]
}else {

############################################################################
#Create matrix for entire population
#Age-Counts from 2020 Census (https://www.census.gov/library/visualizations/interactive/exploring-age-groups-in-the-2020-census.html)
#No immunocompromised
#Unvaccinated prevalence is from CDPH vaccination board (https://covid19.ca.gov/vaccination-progress-data/#overview)
#For vaccinated prevalence, keeping same proportions of 3-dose vs. 4-dose (https://www.cdc.gov/mmwr/volumes/72/wr/mm7207a5.htm#T1_down)

set.seed(488)
age_0_17_cal <- data.frame(individual = c(1:(1103000)),
                           age_group = '0-17 years',
                           num_doses = sample(c('unvax', '3-dose'), 1103000, prob = c(0.574, 0.426), replace = TRUE),
                           prior_inf = rbinom(1103000, 1, 0.919), #updated with pediatric commercial lab seroprevalence (Oct-Dec 2022)
                           immunocompromised = sample(c(0,1,2), 1103000, prob = c(1, 0, 0), replace = TRUE))
set.seed(488)
age_18_49_cal <- data.frame(individual = c((1103000 + 1):(3205000)),
                            age_group = '18-49 years',
                            num_doses = sample(c('unvax', '3-dose'), 2102000, prob = c(0.221, 0.779), replace = TRUE),
                            prior_inf = rbinom(2102000, 1, 0.86), #updated with blood donor seroprevalence (age-adjusted; Oct-Dec 2022)
                            immunocompromised = sample(c(0,1,2), 2102000, prob = c(1, 0, 0), replace = TRUE))

set.seed(488)
age_50_64_cal <- data.frame(individual = c((3205000 + 1):(4158500)),
                            age_group = '50-64 years',
                            num_doses = sample(c('unvax','3-dose', '4-dose'), 953500, prob = c(0.122, 0.6*0.878, 0.4*0.878), replace = TRUE),
                            prior_inf = rbinom(953500, 1, 0.775), #updated with blood donor seroprevalence (Oct-Dec 2022)
                            immunocompromised = sample(c(0,1,2), 953500, prob = c(1, 0, 0), replace = TRUE))

set.seed(488)
age_65_74_cal <- data.frame(individual = c((4158500 + 1):(4531500)),
                            age_group = '65-74 years',
                            num_doses = sample(c('unvax','3-dose', '4-dose'), 373000, prob = c(0.067, 0.6*0.933, 0.4*0.933), replace = TRUE),
                            prior_inf = rbinom(373000, 1, 0.565), #updated with blood donor seroprevalence (Oct-Dec 2022)
                            immunocompromised = sample(c(0,1,2), 373000, prob = c(1, 0, 0), replace = TRUE))

set.seed(488)
age_75_plus_cal <- data.frame(individual = c((4531500  + 1):(5000000)),
                              age_group = '75+ years',
                              num_doses = sample(c('unvax','3-dose', '4-dose'), 468500, prob = c(0.067, 0.6*0.933, 0.4*0.933), replace = TRUE),
                              prior_inf = rbinom(468500, 1, 0.565), #updated with blood donor seroprevalence (Oct-Dec 2022)
                              immunocompromised = sample(c(0,1,2), 468500, prob = c(1, 0, 0), replace = TRUE))

entire_population <- rbind(age_0_17_cal, age_18_49_cal, age_50_64_cal, age_65_74_cal, age_75_plus_cal)

#Clear age-specific dfs
rm(age_0_17_cal, age_18_49_cal, age_50_64_cal, age_65_74_cal, age_75_plus_cal)
############################################################################
#Functions for assigning a 'time-since-last' to each individual

#Calculate time since last vaccine dose or infection
add.weeks= function(date,n) {seq(date, by = paste (n, "weeks"), length = 2)[2]}

time_since_last <- function(df) {
  #Calculate time since last dose and time since last infection
  set.seed(488)
  last_dose_and_inf <- df %>% mutate(time_since_last_dose = ifelse(num_doses == "3-dose", 
                                                                   sample(three_doses_by_week$week,
                                                                          size = sum(num_doses == '3-dose'),
                                                                          prob = three_doses_by_week$perc_doses, 
                                                                          replace = TRUE),
                                                                   ifelse(num_doses == '4-dose',
                                                                          sample(four_doses_by_week$week, 
                                                                                 size = sum(num_doses == '4-dose'),
                                                                                 prob = four_doses_by_week$perc_doses, 
                                                                                 replace = TRUE),
                                                                          '2020-09-01')),
                                     time_since_last_inf = ifelse(prior_inf == 1,
                                                                  sample(as.character(cases_by_week$week), 
                                                                         size = sum(prior_inf == 1),
                                                                         prob = cases_by_week$perc_cases,
                                                                         replace = TRUE),
                                                                  NA))
  #Reinfection for 10% of infected individuals
  set.seed(488)
  reinfection <- last_dose_and_inf %>% filter(prior_inf == 1) %>%  
    sample_frac(.1) %>% 
    mutate(reinf_period = interval(as.Date(time_since_last_inf), as.Date("2022-09-01"),) %/% weeks(1)) %>%
    rowwise() %>% mutate(time_since_last_reinf = ifelse(reinf_period > 0,
                                                        sample(as.character(cases_by_week$week[as.Date(cases_by_week$week) >= add.weeks(as.Date(time_since_last_inf), 3)]),
                                                               size = 1,
                                                               prob = cases_by_week$perc_cases[as.Date(cases_by_week$week) >= add.weeks(as.Date(time_since_last_inf), 3)],
                                                               replace = TRUE),
                                                        NA)) %>% dplyr::select(individual,time_since_last_reinf)                                                          
  
  
  #Merge reinfection data to main df
  merged <- as.data.frame(merge(setDT(last_dose_and_inf), setDT(reinfection), by = c("individual"), all.x = TRUE)) %>%
    mutate(time_since_last_dose_inf = pmax(as.Date(time_since_last_dose), as.Date(time_since_last_inf), as.Date(time_since_last_reinf), as.Date("2020-12-31"), na.rm =  TRUE),
           weeks_since_last_dose_inf = as.numeric(as.character(as.factor(interval(time_since_last_dose_inf,as.Date('2022-12-31')) %/% weeks(1)))),
           weeks_since_last_inf = as.numeric(as.character(as.factor(interval(time_since_last_inf,as.Date('2022-12-31')) %/% weeks(1)))),
           weeks_since_last_reinf = as.numeric(as.character(as.factor(interval(time_since_last_reinf,as.Date('2022-12-31')) %/% weeks(1)))),
           weeks_since_last_dose = as.numeric(as.character(as.factor(interval(time_since_last_dose, as.Date('2022-12-31')) %/% weeks(1))))) %>%
    dplyr::select(-c("time_since_last_dose", "time_since_last_inf", "time_since_last_reinf")) %>%
    rowwise() %>%
    mutate(weeks_since_last_inf = min(weeks_since_last_inf, 104),
           weeks_since_last_reinf = min(weeks_since_last_reinf, 104))

  return(merged)
}

#Run function on population df
entire_pop <- time_since_last(entire_population)

write.csv(entire_pop, "data/clean-data/entire_pop_5mil.csv")

rm(entire_population)
}

############################################################################################
