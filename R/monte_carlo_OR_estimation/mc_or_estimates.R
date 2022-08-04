# final monte carlo estimates of the odds ratios in primary simulation settings
library(tidyverse)

source(here::here("R", "tp_res_functions.R"))
# primary setting 3 -------------------------------------------------------
one_mil <- read_rds(here::here("R", 
                               "monte_carlo_OR_estimation",
                               "simulation_results",
                               "summaryodds_3diff_1milv2.rds"))


onemil_delta <- delta_calc(one_mil$mean_prob_inf[one_mil$current.in == "A"],
                           one_mil$mean_prob_inf[one_mil$current.in == "B"],
                           one_mil$sd_prob_inf[one_mil$current.in == "A"],
                           one_mil$sd_prob_inf[one_mil$current.in == "B"],
                           one_mil$sim_total[one_mil$current.in == "A"])


onemil_delta[1] - 2*onemil_delta[2]
onemil_delta[1] + 2*onemil_delta[2]




# primary setting 0.3 -----------------------------------------------------


o3_diff <- read_rds(here::here("R", 
                               "monte_carlo_OR_estimation",
                               "simulation_results", 
                               "summaryodds_0.3diff.rds"))


o3_delta <- delta_calc(o3_diff$mean_prob_inf[o3_diff$current.in == "A"],
                       o3_diff$mean_prob_inf[o3_diff$current.in == "B"],
                       o3_diff$sd_prob_inf[o3_diff$current.in == "A"],
                       o3_diff$sd_prob_inf[o3_diff$current.in == "B"],
                       o3_diff$sim_total[o3_diff$current.in == "A"])


o3_delta[1] - 2*o3_delta[2]
o3_delta[1] + 2*o3_delta[2]



# primary setting four sevenths -------------------------------------------


o6_diff <- read_rds(here::here("R", 
                               "monte_carlo_OR_estimation",
                               "simulation_results", 
                               "summaryodds_0.6diff.rds"))


o6_delta <- delta_calc(o6_diff$mean_prob_inf[o6_diff$current.in == "A"],
                       o6_diff$mean_prob_inf[o6_diff$current.in == "B"],
                       o6_diff$sd_prob_inf[o6_diff$current.in == "A"],
                       o6_diff$sd_prob_inf[o6_diff$current.in == "B"],
                       o6_diff$sim_total[o6_diff$current.in == "A"])

#solid on this as well

o6_delta[1] - 2*o6_delta[2]
o6_delta[1] + 2*o6_delta[2]


# primary setting 1.75 ----------------------------------------------------

onesevenfive_diff <- read_rds(here::here("R", 
                                         "monte_carlo_OR_estimation",
                                         "simulation_results", 
                                         "summaryodds_1.75diff_1mil.rds"))


onesevenfive_delta <- delta_calc(onesevenfive_diff$mean_prob_inf[onesevenfive_diff$current.in == "A"],
                                 onesevenfive_diff$mean_prob_inf[onesevenfive_diff$current.in == "B"],
                                 onesevenfive_diff$sd_prob_inf[onesevenfive_diff$current.in == "A"],
                                 onesevenfive_diff$sd_prob_inf[onesevenfive_diff$current.in == "B"],
                                 onesevenfive_diff$sim_total[onesevenfive_diff$current.in == "B"])


onesevenfive_delta[1] - 2*onesevenfive_delta[2]
onesevenfive_delta[1] + 2*onesevenfive_delta[2]

