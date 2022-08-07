# running ME2 on primary 1.75

library(tidyverse)
library(rstan)
library(tidybayes)
set.seed(1234)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("tp_res_functions.R")

args <- commandArgs(trailingOnly=TRUE)
print(args)
seed <- as.integer(args[1])
print(seed)

set.seed(as.numeric(seed))

# read in the results -----------------------------------------------------
res_name_suffix <- "16active_8yearsim_1.75diff_res_seed_"

# list all the setting descriptors
setting <- "16active_8year_1.75diff"





file_list <- list.files(pattern = res_name_suffix)


res_list <- map(file_list,  ~read_rds(.x))

prob_inf <- map(res_list, pluck, 3) %>%
  map(~.x %>%
        mutate(numeric_hiv = as.numeric(current.in == "A")))

prob_inf_single <- prob_inf[[seed]]
model_objects_test <- list(N = dim(prob_inf_single)[1],
                                         z = prob_inf_single$tp_infector,
                                         x = prob_inf_single$numeric_hiv,
                                         specificity = .97,
                                         sensitivity = .28)




test_fit <- stan(file ="misclass_logistic_regression.stan",
                                          data = model_objects_test,
                                          seed = 45,
                                          iter = 2000,
                                          chains = 4
)

file_name = file_list[[seed]]

seed_number <- as.numeric(tail(str_extract_all(file_name, pattern = "\\d+")[[1]], 1))
draws <- test_fit %>%
  tidy_draws() %>%
  bind_rows(.id = "sim") %>%
  dplyr::select(sim, beta0, beta1) %>%
  group_by(sim) %>%
  mean_qi() %>%
  mutate(seed = seed_number)

write_csv(draws, str_c("8year_1.75diff_standraws_seed_", seed, ".csv", ""))

