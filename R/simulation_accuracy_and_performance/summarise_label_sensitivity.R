# Summarise re-running of simulation results with different IS cutoffs
# Originally run in separate files
# Simulation groupings are preserved to match results files in simulation_results


# cutoff 0.4 -------------------------------------------------------------

library(tidyverse)
library(rstan)
set.seed(1234)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#source(here::here("code", "R Code", "tp_res_functions.R"))
source("tp_res_functions.R")


# read in the results -----------------------------------------------------
res_name_suffix <- "primary_1.75_res_seed_"

# list all the setting descriptors
setting <- "16active_1.75diff_0.4"





file_list <- list.files(pattern = res_name_suffix)
total_files <- length(file_list)
print(total_files)
file_list <- file_list[sample(1:total_files, 100)]

res_list <- map(file_list,  ~read_rds(.x))




trials <- length(res_list)


# do 0.4 analysis
acc <- tp_acc(res_list, trials, cutoff = 0.4)%>%
  mutate(setting = setting)



specificity <- acc$mean_percent[acc$true_infector == FALSE]
sensitivity <- acc$mean_percent[acc$true_infector == TRUE]


# calculate true values ---------------------------------------------------


true_val <- 1.93
test_summary <- vector(mode = "list", length = 1)
for (i in 1:length(res_list)) {
  res <- res_list
  num <- trials
  spec <- specificity
  sens <- sensitivity
  test_summary[[i]] <- final_summary(res,
                                     spec,
                                     sens,
                                     num,
                                     cluster = TRUE,
                                     true_val,
                                     stan = FALSE,
                                     cutoff = 0.4) %>%
    mutate(setting = setting,
           true_val = true_val)
  
}

# ME2 run separately for computational reasons
# See simulation_results/primary_1.75_0.4sens_stan_template.R
res_name_suffix <- "8year_1.75diff_0.4sens_standraws_seed_"



file_list <- list.files(pattern = res_name_suffix)

draw_list <- map(file_list, ~read_csv(.x))

draws <- draw_list %>%
  bind_rows(.id = "sim")


misclass_results <- draws %>%
  mutate(OR = exp(beta1),
         OR_Low = exp(beta1.lower),
         OR_High = exp(beta1.upper))


stan_res_summary <- misclass_results %>%
  mutate(Coverage = true_val >= OR_Low & true_val <= OR_High,
         reject_low = OR_Low > 1,
         reject_high = OR_High < 1,
         reject = !(OR_Low < 1 & 1 < OR_High),
         bias = OR - true_val,
         sq_error = (OR - true_val)^2,
         CI_width = OR_High - OR_Low,
         diff = abs(OR - true_val)) %>%
  summarise(
    mean_OR = mean(OR),
    mean_lb = mean(OR_Low),
    mean_ub = mean(OR_High),
    mean_diff = mean(diff),
    mean_bias = mean(bias),
    sq_MSE = sqrt(mean(sq_error)),
    percent_reject = mean(reject),
    mean_ci_width = mean(CI_width),
    percent_reject_low = mean(reject_low),
    percent_reject_high = mean(reject_high),
    percent_coverage = mean(Coverage),
    type = "Stan",
    mean_num_samples = unique(test_summary[[1]]$mean_num_samples)
  ) %>%
  mutate(num_trials = trials,
         setting = setting,
         true_val = true_val)

final_summary <- rbind(test_summary[[1]], stan_res_summary)

write_rds(acc, "primary_1.75_acc_tablesv2_0.4sens.rds")
write_rds(final_summary, "primary_1.75_or_tables_0.4sens.rds")




# cutoff 0.5 --------------------------------------------------------------

setting <- "16active_1.75diff_0.5"


acc <- tp_acc(res_list, trials, cutoff = 0.5)%>%
  mutate(setting = setting)



specificity <- acc$mean_percent[acc$true_infector == FALSE]
sensitivity <- acc$mean_percent[acc$true_infector == TRUE]


# calculate true values ---------------------------------------------------


true_val <- 1.93
test_summary <- vector(mode = "list", length = 1)
for (i in 1:length(res_list)) {
  res <- res_list
  num <- trials
  spec <- specificity
  sens <- sensitivity
  test_summary[[i]] <- final_summary(res,
                                     spec,
                                     sens,
                                     num,
                                     cluster = TRUE,
                                     true_val,
                                     stan = FALSE,
                                     cutoff = 0.5) %>%
    mutate(setting = setting,
           true_val = true_val)
  
}

# ME2 run separately for computational efficiency
# See simulation_results/primary_1.75_0.5sens_stan_template.R
res_name_suffix <- "8year_1.75diff_0.5sens_standraws_seed_"



file_list <- list.files(pattern = res_name_suffix)

draw_list <- map(file_list, ~read_csv(.x))

draws <- draw_list %>%
  bind_rows(.id = "sim")


misclass_results <- draws %>%
  mutate(OR = exp(beta1),
         OR_Low = exp(beta1.lower),
         OR_High = exp(beta1.upper))


stan_res_summary <- misclass_results %>%
  mutate(Coverage = true_val >= OR_Low & true_val <= OR_High,
         reject_low = OR_Low > 1,
         reject_high = OR_High < 1,
         reject = !(OR_Low < 1 & 1 < OR_High),
         bias = OR - true_val,
         sq_error = (OR - true_val)^2,
         CI_width = OR_High - OR_Low,
         diff = abs(OR - true_val)) %>%
  summarise(
    mean_OR = mean(OR),
    mean_lb = mean(OR_Low),
    mean_ub = mean(OR_High),
    mean_diff = mean(diff),
    mean_bias = mean(bias),
    sq_MSE = sqrt(mean(sq_error)),
    percent_reject = mean(reject),
    mean_ci_width = mean(CI_width),
    percent_reject_low = mean(reject_low),
    percent_reject_high = mean(reject_high),
    percent_coverage = mean(Coverage),
    type = "Stan",
    mean_num_samples = unique(test_summary[[1]]$mean_num_samples)
  ) %>%
  mutate(num_trials = trials,
         setting = setting,
         true_val = true_val)

final_summary <- rbind(test_summary[[1]], stan_res_summary)

write_rds(acc, "primary_1.75_acc_tablesv2_0.5sens.rds")
write_rds(final_summary, "primary_1.75_or_tables_0.5sens.rds")


# cutoff 0.7 --------------------------------------------------------------


setting <- "16active_1.75diff_0.7"
acc <- tp_acc(res_list, trials, cutoff = 0.7)%>%
  mutate(setting = setting)



specificity <- acc$mean_percent[acc$true_infector == FALSE]
sensitivity <- acc$mean_percent[acc$true_infector == TRUE]


# calculate true values ---------------------------------------------------


true_val <- 1.93
test_summary <- vector(mode = "list", length = 1)
for (i in 1:length(res_list)) {
  res <- res_list
  num <- trials
  spec <- specificity
  sens <- sensitivity
  test_summary[[i]] <- final_summary(res,
                                     spec,
                                     sens,
                                     num,
                                     cluster = TRUE,
                                     true_val,
                                     stan = FALSE,
                                     cutoff = 0.7) %>%
    mutate(setting = setting,
           true_val = true_val)
  
}

# ME2 run separately for computational efficiency
# see simulation_results/primary_1.75_0.7sens_stan_template.R
res_name_suffix <- "8year_1.75diff_0.7sens_standraws_seed_"


file_list <- list.files(pattern = res_name_suffix)

draw_list <- map(file_list, ~read_csv(.x))

draws <- draw_list %>%
  bind_rows(.id = "sim")


misclass_results <- draws %>%
  mutate(OR = exp(beta1),
         OR_Low = exp(beta1.lower),
         OR_High = exp(beta1.upper))


stan_res_summary <- misclass_results %>%
  mutate(Coverage = true_val >= OR_Low & true_val <= OR_High,
         reject_low = OR_Low > 1,
         reject_high = OR_High < 1,
         reject = !(OR_Low < 1 & 1 < OR_High),
         bias = OR - true_val,
         sq_error = (OR - true_val)^2,
         CI_width = OR_High - OR_Low,
         diff = abs(OR - true_val)) %>%
  summarise(
    mean_OR = mean(OR),
    mean_lb = mean(OR_Low),
    mean_ub = mean(OR_High),
    mean_diff = mean(diff),
    mean_bias = mean(bias),
    sq_MSE = sqrt(mean(sq_error)),
    percent_reject = mean(reject),
    mean_ci_width = mean(CI_width),
    percent_reject_low = mean(reject_low),
    percent_reject_high = mean(reject_high),
    percent_coverage = mean(Coverage),
    type = "Stan",
    mean_num_samples = unique(test_summary[[1]]$mean_num_samples)
  ) %>%
  mutate(num_trials = trials,
         setting = setting,
         true_val = true_val)

final_summary <- rbind(test_summary[[1]], stan_res_summary)

write_rds(acc, "primary_1.75_acc_tablesv2_0.7sens.rds")
write_rds(final_summary, "primary_1.75_or_tables_0.7sens.rds")


