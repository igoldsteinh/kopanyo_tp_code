# summarise secondary quadruple sampling density simulation results
library(tidyverse)
library(rstan)
set.seed(1234)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("tp_res_functions.R")


# read in the results -----------------------------------------------------
res_name_suffix <- list("secondary_quadruple_density_res_seed_")

setting <- list("quadruple_density")

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

 
 glm_tables <- map(res_list, glm_tables, 1.94) %>%
   map2(setting, ~.x %>% mutate(setting = .y))
 
# calculate true values ---------------------------------------------------
# 
# true_val <- c(1.94)
# test_summary <- vector(mode = "list", length = length(res_list))
# for (i in 1:length(res_list)) {
#   res <- res_list[[i]]
#   num <- trials[i]
#   spec <- specificity[[i]]
#   sens <- sensitivity[[i]]
#   test_summary[[i]] <- final_summary(res,
#                                      spec,
#                                      sens,
#                                      num,
#                                      cluster = TRUE,
#                                      true_val[i]) %>%
#                        mutate(setting = setting[[i]],
#                               true_val = true_val[i])
# 
# }
# 
# write_rds(acc, "quadruple_density_acc_tables.rds")
# write_rds(test_summary, "quadruple_density_or_tables.rds")
#write_rds(trials, "trials.rds")
 
 write_rds(glm_tables, "quadruple_density_glmtables.rds")