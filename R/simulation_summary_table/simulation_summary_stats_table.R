# File used to create values reported in 
# the table "Simulation Summary Statistics for 8 Year Simulation."
# Simulations run in primary_1.75_fullsims.R
library(tidyverse)
library(rstan)
set.seed(1234)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#source(here::here("code", "R Code", "tp_res_functions.R"))
source(here::here("R", "tp_res_functions.R"))

res_name_suffix <- "16active_8yearsim_1.75diff_fullsim_seed_"

# res_name_suffix <- list("fulloutbreak_sim_res_seed_",
#                         "fulloutbreak_highsample_sim_res_seed_")
# list all the setting descriptors
setting <- "1.75diff"

# address <- "C:/Users/fiddl/Documents/kopanyo-archived-phylo/code/R Code/fulloutbreak_sim"

# file_list <- map(res_name_suffix, ~list.files(address, 
#                         pattern = .x))

file_list <- list.files(pattern = res_name_suffix)

print(file_list)
exclude <- c("16active_8yearsim_1.75diff_fullsim_seed_100.rds",
               "16active_8yearsim_1.75diff_fullsim_seed_22.rds",
               "16active_8yearsim_1.75diff_fullsim_seed_37.rds",
               "16active_8yearsim_1.75diff_fullsim_seed_5.rds",
               "16active_8yearsim_1.75diff_fullsim_seed_52.rds",
               "16active_8yearsim_1.75diff_fullsim_seed_61.rds",
               "16active_8yearsim_1.75diff_fullsim_seed_64.rds",
               "16active_8yearsim_1.75diff_fullsim_seed_67.rds",
               "16active_8yearsim_1.75diff_fullsim_seed_85.rds")

file_list <- file_list[-which(file_list %in% exclude)]
total_files <- length(file_list)
file_list <- file_list[sample(1:total_files, 100)]
print(file_list)
print(length(file_list))
# # res_list <- map(file_list, ~map(.x, ~read_rds(here::here("code", 
#                                             "R Code", 
#                                             "fulloutbreak_sim", 
#                                             .x))))


res_list <- map(file_list, ~read_rds(.x))


probinf_list <- map(res_list, pluck, "samples")




# paper table -------------------------------------------------------------
outputs_list <- map(res_list, pluck, "outputs")
outputs_list <- map(outputs_list, ~map(.x, ~.x %>%
                                         mutate(true_infector = hosts.ID %in% inf.by)))

outputs_frames <- map(outputs_list, ~bind_rows(., .id = "sim"))


# do al calc real quick
probs_frames <- map(outputs_frames, ~.x %>% filter(active_time > 0))

all_probs <- calc_all_probs(probs_frames) %>%
  mutate(setting = setting)

big_frame <- outputs_frames %>% bind_rows(.id = "simsim")

big_frame$active_time[big_frame$active_time < 0] <- 0

cluster_size <- big_frame %>%
  ungroup() %>%
  group_by(simsim, sim ) %>%
  summarise(cluster_size = n()) %>%
  ungroup() %>%
  summarise(quantiles = quantile(cluster_size, probs = c(0.025, .5, .975)))


summary_stats_by_sim <- big_frame %>%
  ungroup() %>%
  filter(active_time > 0) %>%
  mutate(incub_length = difftime - active_time) %>%
  summarise(quantiles_incub = quantile(incub_length, probs = c(0.025, .5, .975)),
            quantiles_active = quantile(active_time, probs = c(0.025, .5, .975)))


# count number sampled from each cluster
prob_inf <- map(res_list, pluck, "samples") %>%
  map(~.x %>%
        mutate(numeric_hiv = as.numeric(current.in == "A")))

prob_inf_frame <- prob_inf %>% 
  bind_rows(.id = "simsim")

count_total <- prob_inf_frame %>%
  group_by(simsim, sim) %>%
  summarise(num_sampled = n()) %>%
  ungroup() %>%
  summarise(
    quantiles = quantile(num_sampled, probs = c(0.025, .5, .975)))

# mean tree heights -------------------------------------------------------

library(phytools)
tree_heights <- purrr::map(res_list, pluck, "sample_tree") %>%
  purrr::map(~purrr::map(.x, nodeHeights)) %>%
  purrr::map(~purrr::map(.x, max)) %>%
  purrr::map(unlist) %>%
  unlist() 



mean(tree_heights)
sd(tree_heights)

quantile(tree_heights, c(0.025, 0.5, 0.975))
write_csv(cluster_size, "1.75diff_clustersize.csv")
write_csv(summary_stats_by_sim, "1.75diff_summarytimes.csv")
write_csv(count_total, "1.75diff_samplecount.csv")
write_rds(all_probs, "1.75diff_allprobs.csv")


