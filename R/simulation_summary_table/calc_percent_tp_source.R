# file for calculating what percent of sampled cases
# are identified as infection sources using TP

# Summarise primary and secondary simulations
# Originally run in separate files
# Simulation groupings are preserved to match results files in simulation_results

library(tidyverse)
library(rstan)
set.seed(1234)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#source(here::here("code", "R Code", "tp_res_functions.R"))
source(here::here("tp_res_functions.R"))




# primary settings 1.75 and 4/7 -------------------------------------------


res_name_suffix <- list("16active_8yearsim_1.75diff_res_seed_")

# list all the setting descriptors
setting <- list("primary_1.75")

file_list <- map(res_name_suffix, ~list.files(pattern = .x))


total_files <- map(file_list, ~length(.x))
print(total_files)
file_list <- map2(file_list, total_files, ~.x[sample(1:.y, 100)])

res_list <- map(file_list, ~map(.x, ~read_rds(.x)))

res <- res_list[[1]]

prob_inf_frame <- map(res, pluck, "final_inf_frame")
                  
  
percent_tp_source <- map(prob_inf_frame, ~.x %>%
                                         ungroup() %>%
                                         summarise(percent_tp = sum(tp_infector)/n())) %>%
                     bind_rows() %>%
                     pull()



print(median(percent_tp_source))
print(quantile(percent_tp_source, c(0.025, 0.975)))
# test <- read_rds("R/simulation_summary_table/16active_8yearsim_1.75diff_res_seed_105.rds")
# 
# test_frame <- pluck(test, "final_inf_frame")

