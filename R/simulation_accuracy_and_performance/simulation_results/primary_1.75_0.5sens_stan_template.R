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
res_name_suffix <- "16active_8yearsim_1.75diff_res_seed_"

# res_name_suffix <- list("fulloutbreak_sim_res_seed_",
#                         "fulloutbreak_highsample_sim_res_seed_")
# list all the setting descriptors
setting <- "16active_1.75diff_0.4"



# address <- "C:/Users/fiddl/Documents/kopanyo-archived-phylo/code/R Code/fulloutbreak_sim"


file_list <- list.files(pattern = res_name_suffix)
# file_list <- list.files(path = address, pattern = res_name_suffix)

# res_list <- map(file_list, ~read_rds(here::here("code",
#                                            "R Code",
#                                            "fulloutbreak_sim",
#                                            .x)))

total_files <- length(file_list)
print(total_files)
file_list <- file_list[sample(1:total_files, 100)]

res_list <- map(file_list,  ~read_rds(.x))




trials <- length(res_list)


# do 0.4 analysis
acc <- tp_acc(res_list, trials, cutoff = 0.5)%>%
  mutate(setting = setting)



specificity <- acc$mean_percent[acc$true_infector == FALSE]
sensitivity <- acc$mean_percent[acc$true_infector == TRUE]
# test_summary <- vector(mode = "list", length = length(res_list))

# test_summary <- vector(mode = "list", length = length(res_list))
prob_inf <- map(res_list, pluck, 3) %>%
  map(~.x %>%
        mutate(numeric_hiv = as.numeric(current.in == "A"),
               tp_infector = prob_inf > 0.5))

prob_inf_single <- prob_inf[[seed]]
model_objects_test <- list(N = dim(prob_inf_single)[1],
                                         z = prob_inf_single$tp_infector,
                                         x = prob_inf_single$numeric_hiv,
                                         specificity = specificity,
                                         sensitivity = sensitivity)




test_fit <- stan(file ="misclass_logistic_regression.stan",
                                          data = model_objects_test,
                                          seed = 45,
                                          iter = 2000,
                                          chains = 4
)
draws <- test_fit %>%
  tidy_draws() %>%
  bind_rows(.id = "sim") %>%
  dplyr::select(sim, beta0, beta1) %>%
  group_by(sim) %>%
  mean_qi()

write_csv(draws, str_c("8year_1.75diff_0.5sens_standraws_seed_", seed, ".csv", ""))

