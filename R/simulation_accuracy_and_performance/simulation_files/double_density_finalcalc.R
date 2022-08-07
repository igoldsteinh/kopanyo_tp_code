# summarise secondary increase density simulation

library(tidyverse)
library(rstan)
set.seed(1234)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("tp_res_functions.R")


# read in the results -----------------------------------------------------
res_name_suffix <- "32active_8yearsim_1.75diff_res_seed_"

setting <- "double_density"
                
              



file_list <- list.files(pattern = res_name_suffix)

total_files <- length(file_list)
print(total_files)
file_list <- file_list[sample(1:total_files, 100)]

seed_list <- map(file_list, ~as.numeric(tail(str_extract_all(.x, pattern = "\\d+")[[1]], 1))) %>%
            unlist()
print(seed_list)
res_list <- map(file_list,  ~read_rds(.x))

trials <- length(res_list)
acc <- tp_acc(res_list, trials)%>%
        mutate(setting = setting)



specificity <- acc$mean_percent[acc$true_infector == FALSE]
sensitivity <- acc$mean_percent[acc$true_infector == TRUE]


# calculate true values ---------------------------------------------------

true_val <- 1.94
test_summary <- vector(mode = "list", length = 1)
for (i in 1:length(res_name_suffix)) {
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


res_name_suffix <- "double_density_standraws_seed_"



file_list <- list.files(pattern = res_name_suffix)

draw_list <- map(file_list, ~read_csv(.x))

draws <- draw_list %>%
  bind_rows(.id = "sim") %>%
  filter(seed %in% seed_list)

length(unique(draws$seed))

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
write_rds(acc, "double_density_acc_tables.rds")
write_rds(final_summary, "double_density_sim_or_tables.rds")
