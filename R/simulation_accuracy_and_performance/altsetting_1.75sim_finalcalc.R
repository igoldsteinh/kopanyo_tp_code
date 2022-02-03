# this file summarises pipeline performance
# across many simulation settings
library(tidyverse)
library(rstan)
set.seed(1234)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
 #source(here::here("code", "R Code", "tp_res_functions.R"))
source("tp_res_functions.R")


# read in the results -----------------------------------------------------
res_name_suffix <- list("32active_8yearsim_1.75diff_res_seed_")

# res_name_suffix <- list("fulloutbreak_sim_res_seed_",
#                         "fulloutbreak_highsample_sim_res_seed_")
# list all the setting descriptors
setting <- list("32active_8yearsim_1.75"
              )

# address <- "C:/Users/fiddl/Documents/kopanyo-archived-phylo/code/R Code/fulloutbreak_sim"

# file_list <- map(res_name_suffix, ~list.files(address, 
#                         pattern = .x))

file_list <- map(res_name_suffix, ~list.files(pattern = .x))
                 
# # res_list <- map(file_list, ~map(.x, ~read_rds(here::here("code", 
#                                             "R Code", 
#                                             "fulloutbreak_sim", 
#                                             .x))))

total_files <- map(file_list, ~length(.x))
print(total_files)
file_list <- map2(file_list, total_files, ~.x[sample(1:.y, 100)])

 res_list <- map(file_list, ~map(.x, ~read_rds(.x)))

trials <- map(res_list, length)
acc <- map2(res_list, trials, ~tp_acc(.x, .y)) %>%
       map2(setting, ~.x %>% mutate(setting = .y))

specificity <- map(acc, ~.x$mean_percent[.x$true_infector==FALSE])
sensitivity <- map(acc, ~.x$mean_percent[.x$true_infector==TRUE])
# test_summary <- vector(mode = "list", length = length(res_list))

glm_tables <- map(res_list, glm_tables, 1.93) %>%
              map2(setting, ~.x %>% mutate(setting = .y))

prob_inf_list <- map(res_list, ~.x %>% map(pluck, "final_inf_frame"))

count_labels <- map(prob_inf_list, count_2x2) %>%
  map2(setting, ~.x %>% mutate(setting = .y))
# calculate true values ---------------------------------------------------


true_val <- c(1.93,1.93, 1.93)
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

write_rds(acc, "altsetting_1.75diff_acc_tablesv1.rds")
# write_rds(glm_tables, "altsetting_1.75diff_glmtables.rds")
# write_rds(count_labels, "altsetting_1.75diff_2x2tables.rds")

write_rds(test_summary, "altsetting_1.75diff_sim_or_tablesv1.rds")
