# summarising odds
library(tidyverse)
#source(here::here("code", "R Code", "tp_res_functions.R"))
source("tp_res_functions.R")
res_name_suffix <- "16active_8yearsim_0.6diff_odds_seed_"

file_list <- list.files(pattern = res_name_suffix)



 # res_list <- map(file_list, ~map(.x, ~read_rds(here::here("code",
 #                                            "R Code",
 #                                            "odds_estimation_files",
 #                                            .x))))


res_list <- map(file_list,  ~read_rds(.x))


summary_odds <- res_list %>%
  bind_rows(.id = "sim_group") %>%
  group_by(sim_group, sim) %>%
  mutate(num_rows = n()) %>%
  filter(num_rows == 2) %>%
  group_by(current.in)%>% 
  summarise(mean_prob_inf = mean(prob_inf),
            sim_total = n(),
            mean_total = mean(num_total),
            sd_prob_inf = sd(prob_inf),
            se_prob_inf = sd(prob_inf)/sqrt(sim_total)
  )



write_rds(summary_odds, "summaryodds_0.6diff.rds")