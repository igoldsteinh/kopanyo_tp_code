# running stan on 3diff

library(tidyverse)
library(rstan)
library(tidybayes)
set.seed(1234)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#source(here::here("code", "R Code", "tp_res_functions.R"))
source("tp_res_functions.R")

args <- commandArgs(trailingOnly=TRUE)
print(args)
seed <- as.integer(args[1])
print(seed)

set.seed(as.numeric(seed))

# read in the results -----------------------------------------------------
res_name_suffix <- "16active_8yearsim_samediff_res_seed_"

# res_name_suffix <- list("fulloutbreak_sim_res_seed_",
#                         "fulloutbreak_highsample_sim_res_seed_")
# list all the setting descriptors
setting <- "16active_8year_samediff"



# address <- "C:/Users/fiddl/Documents/kopanyo-archived-phylo/code/R Code/fulloutbreak_sim"


file_list <- list.files(pattern = res_name_suffix)
# file_list <- list.files(path = address, pattern = res_name_suffix)

# res_list <- map(file_list, ~read_rds(here::here("code",
#                                            "R Code",
#                                            "fulloutbreak_sim",
#                                            .x)))


res_list <- map(file_list,  ~read_rds(.x))

# test_summary <- vector(mode = "list", length = length(res_list))
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

write_csv(draws, str_c("16active_8year_samediff_standraws_seed_", seed, ".csv", ""))

