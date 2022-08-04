library(tidyverse)

test <- read_rds(here::here("R", "simulation_accuracy_and_performance",
                            "simulation_results", "16active_8yearsim_0.3diff_res_seed_8.rds"))
trees <- test[["sample_tree"]]



library(phytools)
tree_heights <- pluck(test, "sample_tree") %>%
  purrr::map(nodeHeights) %>%
  purrr::map(max) %>%
  purrr::map(unlist) %>%
  unlist() 
