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


res_name_suffix <- list("primary_1.75_res_seed_",
                        "primary_foursevenths_res_seed_")

# list all the setting descriptors
setting <- list("primary_1.75",
                "primary_foursevenths"
)

file_list <- map(res_name_suffix, ~list.files(pattern = .x))


total_files <- map(file_list, ~length(.x))
print(total_files)
file_list <- map2(file_list, total_files, ~.x[sample(1:.y, 100)])

res_list <- map(file_list, ~map(.x, ~read_rds(.x)))

trials <- map(res_list, length)
acc <- map2(res_list, trials, ~tp_acc(.x, .y))%>%
  map2(setting, ~.x %>% mutate(setting = .y))
# 
specificity <- map(acc, ~.x$mean_percent[.x$true_infector==FALSE])
sensitivity <- map(acc, ~.x$mean_percent[.x$true_infector==TRUE])
test_summary <- vector(mode = "list", length = length(res_list))


# calculate results ---------------------------------------------------


true_val <- c(1.94,0.52)
test_summary <- vector(mode = "list", length = length(res_list))
for (i in 1:length(res_list)) {
  res <- res_list[[i]]
  num <- trials[i]
  spec <- specificity[[i]]
  sens <- sensitivity[[i]]
  test_summary[[i]] <- final_summary(res,
                                     spec,
                                     sens,
                                     num,
                                     cluster = TRUE,
                                     true_val[i]) %>%
    mutate(setting = setting[[i]],
           true_val = true_val[i])
  
}

write_rds(acc, "primary_1.75_foursevenths_acc_tablesv1.rds")
write_rds(test_summary, "primary_1.75_foursevenths_or_tablesv1.rds")


rm(res_list)
# primary simulation setting 3 --------------------------------------------


# read in the results -----------------------------------------------------
res_name_suffix <- "primary_3_res_seed_"
setting <- "primary_3"



file_list <- list.files(pattern = res_name_suffix)
total_files <- length(file_list)
print(total_files)
file_list <- file_list[sample(1:total_files, 100)]

res_list <- map(file_list,  ~read_rds(.x))

trials <- length(res_list)
acc <- tp_acc(res_list, trials)%>%
  mutate(setting = setting)



specificity <- acc$mean_percent[acc$true_infector == FALSE]
sensitivity <- acc$mean_percent[acc$true_infector == TRUE]


# calculate results ---------------------------------------------------


true_val <- 3.53
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
                                     stan = FALSE) %>%
    mutate(setting = setting,
           true_val = true_val)
  
}

# ME2 models run separately in primary_3_stan_template.R for computational efficiency
# see simulation_results/primary_3_stan_template.R
res_name_suffix <- "primary_3_standraws_seed_"



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
write_rds(acc, "primary_3_acc_tablesv2.rds")
write_rds(final_summary, "primary_3_or_tablesv2.rds")

rm(res_list)


# primary simulation settings 1 and 0.3 ----------------------------------------------
# read in the results -----------------------------------------------------
res_name_suffix <- list("primary_1_res_seed_",
                        "primary_0.3_res_seed_")

# list all the setting descriptors
setting <- list("primary_1",
                "primary_0.3"
)


file_list <- map(res_name_suffix, ~list.files(pattern = .x))


total_files <- map(file_list, ~length(.x))
print(total_files)
file_list <- map2(file_list, total_files, ~.x[sample(1:.y, 100)])

res_list <- map(file_list, ~map(.x, ~read_rds(.x)))

trials <- map(res_list, length)
acc <- map2(res_list, trials, ~tp_acc(.x, .y))%>%
  map2(setting, ~.x %>% mutate(setting = .y))

specificity <- map(acc, ~.x$mean_percent[.x$true_infector==FALSE])
sensitivity <- map(acc, ~.x$mean_percent[.x$true_infector==TRUE])
test_summary <- vector(mode = "list", length = length(res_list))


# calculate results ---------------------------------------------------

true_val <- c(1, 0.25)
test_summary <- vector(mode = "list", length = length(res_list))
for (i in 1:length(res_list)) {
  res <- res_list[[i]]
  num <- trials[i]
  spec <- specificity[[i]]
  sens <- sensitivity[[i]]
  test_summary[[i]] <- final_summary(res,
                                     spec,
                                     sens,
                                     num,
                                     cluster = TRUE,
                                     true_val[i]) %>%
    mutate(setting = setting[[i]],
           true_val = true_val[i])
  
}

write_rds(acc, "primary_1_0.3_acc_tables.rds")
write_rds(test_summary, "primary_1_0.3_or_tables.rds")


# Secondary setting increase sampling window ------------------------------


# read in the results -----------------------------------------------------
res_name_suffix <- list("secondary_increase_window_res_seed_")

setting <- list("increase window")


file_list <- map(res_name_suffix, ~list.files(pattern = .x))

res_list <- map(file_list, ~map(.x, ~read_rds(.x)))

trials <- map(res_list, length)
acc <- map2(res_list, trials, ~tp_acc(.x, .y)) %>%
  map2(setting, ~.x %>% mutate(setting = .y))

specificity <- map(acc, ~.x$mean_percent[.x$true_infector==FALSE])
sensitivity <- map(acc, ~.x$mean_percent[.x$true_infector==TRUE])

# calculate results ---------------------------------------------------


true_val <- c(1.94,1.94, 1.94)
test_summary <- vector(mode = "list", length = length(res_list))
for (i in 1:length(res_list)) {
  res <- res_list[[i]]
  num <- trials[i]
  spec <- specificity[[i]]
  sens <- sensitivity[[i]]
  test_summary[[i]] <- final_summary(res,
                                     spec,
                                     sens,
                                     num,
                                     cluster = TRUE,
                                     true_val[i],
                                     stan = FALSE) %>%
    mutate(setting = setting[[i]],
           true_val = true_val[i])
  
}

# ME2 models run separately in primary_3_stan_template.R for computational efficiency
# see simulation_results/primary_3_stan_template.R

res_name_suffix <- "secondary_increase_window_standraws_seed_"




file_list <- list.files(pattern = res_name_suffix)

draw_list <- map(file_list, ~read_csv(.x))

draws <- draw_list %>%
  bind_rows(.id = "sim")


misclass_results <- draws %>%
  mutate(OR = exp(beta1),
         OR_Low = exp(beta1.lower),
         OR_High = exp(beta1.upper))


stan_res_summary <- misclass_results %>%
  mutate(Coverage = true_val[1] >= OR_Low & true_val[1] <= OR_High,
         reject_low = OR_Low > 1,
         reject_high = OR_High < 1,
         reject = !(OR_Low < 1 & 1 < OR_High),
         bias = OR - true_val[1], 
         sq_error = (OR - true_val[1])^2,
         CI_width = OR_High - OR_Low,
         diff = abs(OR - true_val[1])) %>%
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
         true_val = true_val[1])

final_summary <- rbind(test_summary[[1]], stan_res_summary)

write_rds(acc, "secondary_increase_window_acc_tablesv1.rds")

write_rds(final_summary, "secondary_increase_window_or_tablesv1.rds")


rm(res_list)
# secondary setting increase sampling density -----------------------------

# read in the results -----------------------------------------------------
res_name_suffix <- list("secondary_increase_density_res_seed_")
setting <- list("increase density"
)

file_list <- map(res_name_suffix, ~list.files(pattern = .x))

total_files <- map(file_list, ~length(.x))
print(total_files)
file_list <- map2(file_list, total_files, ~.x[sample(1:.y, 100)])

res_list <- map(file_list, ~map(.x, ~read_rds(.x)))

trials <- map(res_list, length)
acc <- map2(res_list, trials, ~tp_acc(.x, .y)) %>%
  map2(setting, ~.x %>% mutate(setting = .y))

specificity <- map(acc, ~.x$mean_percent[.x$true_infector==FALSE])
sensitivity <- map(acc, ~.x$mean_percent[.x$true_infector==TRUE])
# calculate true values ---------------------------------------------------


true_val <- c(1.93)
test_summary <- vector(mode = "list", length = length(res_list))
for (i in 1:length(res_list)) {
  res <- res_list[[i]]
  num <- trials[i]
  spec <- specificity[[i]]
  sens <- sensitivity[[i]]
  test_summary[[i]] <- final_summary(res,
                                     spec,
                                     sens,
                                     num,
                                     cluster = TRUE,
                                     true_val[i]) %>%
    mutate(setting = setting[[i]],
           true_val = true_val[i])
  
}

write_rds(acc, "secondary_increase_density_acc_tablesv1.rds")
write_rds(test_summary, "secondary_increase_density_or_tablesv1.rds")

rm(res_list)
# secondary increase clusters setting -------------------------------------

# read in the results -----------------------------------------------------
res_name_suffix <- list("secondary_increase_clusters_res_seed_"
)
setting <- list(
  "increase clusters"
  
)

file_list <- map(res_name_suffix, ~list.files(pattern = .x))


total_files <- map(file_list, ~length(.x))
print(total_files)
file_list <- map2(file_list, total_files, ~.x[sample(1:.y, 100)])

res_list <- map(file_list, ~map(.x, ~read_rds(.x)))

trials <- map(res_list, length)
acc <- map2(res_list, trials, ~tp_acc(.x, .y))

specificity <- map(acc, ~.x$mean_percent[.x$true_infector==FALSE])
sensitivity <- map(acc, ~.x$mean_percent[.x$true_infector==TRUE])
# test_summary <- vector(mode = "list", length = length(res_list))


# calculate results ---------------------------------------------------


true_val <- c(1.94)
test_summary <- vector(mode = "list", length = length(res_list))
for (i in 1:length(res_list)) {
  res <- res_list[[i]]
  num <- trials[i]
  spec <- specificity[[i]]
  sens <- sensitivity[[i]]
  test_summary[[i]] <- final_summary(res, 
                                     spec, 
                                     sens, 
                                     num, 
                                     cluster = TRUE, 
                                     true_val[i]) %>%
    mutate(setting = setting[[i]],
           true_val = true_val[i])
  
}

write_rds(acc, "secondary_increase_clusters_acc_tablesv1.rds")
write_rds(test_summary, "secondary_increase_clusters_or_tablesv1.rds")

