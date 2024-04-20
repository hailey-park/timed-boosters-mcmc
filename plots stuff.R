nonsevere_waning <- melt(read.csv("combined_nonsevere_waning_predictions_weekly_prior_inf_only.csv")[,-1] %>%
                           filter(estimate == "mean",
                                  immunocompromised == 1,
                                  age_group == "65-74 years") %>%
                           rowwise() %>%
                           mutate(plus_20 = max(ve_pred + 0.2, 0),
                                  minus_20 = max(ve_pred - 0.2, 0)) %>%
                           dplyr::select(weeks, ve_pred, plus_20, minus_20), id = "weeks")

nonsevere_waning %>% 
  ggplot(aes(weeks, value * 100, group=variable, color=variable)) + 
  geom_line() + geom_point() +
  ylim(0,100) +
  ggtitle("Nonsevere protection curves\nImmunity Type: Prior Infection Only") +
  ylab("Protective Effectiveness (%)") + 
  xlab("Weeks") +
  labs(color = "Protection") +
  scale_color_discrete(labels = c("Mean", "+20%", "-20%")) 





# adjusted[[2]] %>% 
#   filter(#prior_inf == 1,
#          age_group == "50-64 years",
#          immunocompromised == 0) %>%
#   ggplot(aes(weeks, nonsevere_prior_inf_only_ve_pred)) + 
#   geom_line() + geom_point() +
#   ylim(0,1)
# 
# nonsevere_waning_prior_inf_only %>% 
#   filter(#prior_inf == 1,
#          age_group == "50-64 years",
#          immunocompromised == 0) %>%
#   ggplot(aes(weeks, nonsevere_prior_inf_only_ve_pred)) + 
#   geom_line() + geom_point() +
#   ylim(0,1)

# 
# waning_rate_adjustment_fn <- function(rate_vacc, rate_inf, rate_hybrid, df_hybrid_vacc, df_inf){
#   
#   if(rate_vacc < 0){
#     df_hybrid_vacc <- df_hybrid_vacc %>% 
#       rowwise() %>%
#       mutate(nonsevere_ve_pred = if_else(prior_inf == 0, 
#                                          max(nonsevere_ve_pred - (weeks * rate_vacc/104), 0),
#                                          nonsevere_ve_pred))
#   }
#   
#   if(rate_vacc >= 0) {
#     df_hybrid_vacc <- df_hybrid_vacc %>% 
#       rowwise() %>% 
#       mutate(nonsevere_ve_pred = if_else(prior_inf == 0, 
#                                          max(nonsevere_ve_pred - (weeks * rate_vacc/104), 0),
#                                          nonsevere_ve_pred))
#   }
#   
#   if(rate_inf < 0){
#     df_inf <- df_inf %>% 
#       rowwise() %>%
#       mutate(nonsevere_prior_inf_only_ve_pred = max(nonsevere_prior_inf_only_ve_pred - (weeks * rate_inf/104), 0))
#   }
#   
#   if(rate_inf >= 0) {
#     df_inf <- df_inf %>% 
#       rowwise() %>% 
#       mutate(nonsevere_prior_inf_only_ve_pred = max(nonsevere_prior_inf_only_ve_pred - (weeks * rate_inf/104), 0))
#     
#   }
#   
#   if(rate_hybrid < 0){
#     df_hybrid_vacc <- df_hybrid_vacc %>% 
#       rowwise() %>%
#       mutate(nonsevere_ve_pred = if_else(prior_inf == 1, 
#                                          max(nonsevere_ve_pred - (weeks * rate_hybrid/104), 0),
#                                          nonsevere_ve_pred))
#   }
#   
#   if(rate_hybrid >= 0) {
#     df_hybrid_vacc <- df_hybrid_vacc %>% 
#       rowwise() %>% 
#       mutate(nonsevere_ve_pred = if_else(prior_inf == 1, 
#                                          max(nonsevere_ve_pred - (weeks * rate_hybrid/104), 0),
#                                          nonsevere_ve_pred))
#   }
#   
#   return(list(df_hybrid_vacc, df_inf))
# }
# 
# adjusted_waning <- waning_rate_adjustment_fn(pe_nonsevere_rate_vacc, pe_nonsevere_rate_inf, pe_nonsevere_rate_hybrid, nonsevere_waning, nonsevere_waning_prior_inf_only)
# 



hosp_data <- read.csv("data/raw-data/Weekly_Rates_of_Laboratory-Confirmed_COVID-19_Hospitalizations_from_the_COVID-NET_Surveillance_System_20240325.csv")

weekly_inc_gen_pop <- hosp_data %>% 
  mutate(week = as.Date(str_sub(Week.ending.date, 1, 10), format = "%m/%d/%Y"),
         weeks_since = as.numeric(interval(as.Date('2022-12-01'), week) %/% weeks(1)),
         Age.Category = if_else(Age.Category == "0-17 years (Children)", "0-17 years", Age.Category)) %>%
  filter(State == "COVID-NET",
         Season %in% c("2021-22", "2022-23", "2023-24"),
         Age.Category %in% c("All", "0-17 years", "18-49 years", "50-64 years", "65-74 years", "75+ years"),
         Sex == "All",
         Race == "All",
         week >= as.Date("2022-12-01") & week <= as.Date("2023-12-31")) %>%
  dplyr::select(week, weeks_since, Age.Category, Rate) %>%
  rename(age_group = Age.Category,
         observed_inc = Rate) 

write.csv(weekly_inc_gen_pop, "data/clean-data/observed-inc-estimates-csv")
