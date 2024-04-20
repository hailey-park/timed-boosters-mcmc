###################################################################################################
#Title: Cleaning Case Data and Vaccine Data Used for 'Time Since Last' Estimation
#Author: Hailey Park
#Date: March 27, 2023
###################################################################################################

rm(list=ls())

setwd(here::here())

#Load libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(lubridate)


#Read in data
#Case data from this website (https://data.cdc.gov/Case-Surveillance/Weekly-United-States-COVID-19-Cases-and-Deaths-by-/pwn4-m3yp)
#Vaccine data from this website (https://data.chhs.ca.gov/dataset/vaccine-progress-dashboard/resource/a4b01fad-f75d-4a8c-a88a-1513a07fa4b2)

cases <- read.csv("data/raw-data/Weekly_United_States_COVID-19_Cases_and_Deaths_by_State_-_ARCHIVED.csv")
vaccines <- read.csv("data/raw-data/covid19vaccinesadministeredbydemographics-booster-archived.csv")


#Cases By Week
cases_by_week <- cases %>% group_by(end_date) %>% summarise(num_cases = sum(new_cases)) %>% mutate(week = as.Date(end_date, format="%m/%d/%Y")) %>%
  filter(week >= as.Date('2020-09-01') & week <= as.Date('2022-12-31')) %>% mutate(perc_cases = round(num_cases/sum(num_cases) * 100, 2))

#Booster Doses by Week
#NOTE:    We are assuming any booster doses received after 3/1/22 is a second booster dose. 
#         We are also excluding anyone from receiving a booster dose 2 months before the simulation 
#         start date (1/1/23).
doses_by_week <- vaccines %>% mutate(week = floor_date(as.Date(administered_date), unit = "week"), administered_date = as.Date(administered_date)) %>% 
  filter((administered_date >= '2021-09-22' & administered_date <= '2022-12-31'), demographic_category == "Age Group", demographic_value %in% c("18-49", "50-64", "65+")) %>%
  mutate(total_booster_doses = total_doses - partially_vaccinated - fully_vaccinated) %>% group_by(week) %>% summarise(total_boosted = sum(booster_recip_count), total_booster_doses = sum(total_booster_doses)) %>%
  dplyr::select(week, total_booster_doses)


four_doses_by_week <- doses_by_week %>% filter(week >= as.Date("2022-03-01") & week <= as.Date("2022-10-01")) %>% mutate(perc_doses = round(total_booster_doses/sum(total_booster_doses), 3))
three_doses_by_week <- doses_by_week %>% filter(week <= as.Date("2022-03-01")) %>% mutate(perc_doses = round(total_booster_doses/sum(total_booster_doses), 3))


#write as .csv
write.csv(cases_by_week, "data/clean-data/cases_by_week.csv")
write.csv(four_doses_by_week, "data/clean-data/four_doses_by_week.csv")
write.csv(three_doses_by_week, "data/clean-data/three_doses_by_week.csv")


