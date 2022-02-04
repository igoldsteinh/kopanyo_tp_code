# Produces visualizations of results when cutoff for IS label is adjusted
# Uses results produced in summarise_label_sensitivity.R
library(tidyverse)
library(xtable)
library(patchwork)

acc_twov1 <- read_rds(here::here("R",
                                 "simulation_accuracy_and_performance",
                                 "simulation_results",
                                 "primary_1.75_foursevenths_acc_tablesv1.rds"))
summary_twov1 <- read_rds(here::here("R",
                                     "simulation_accuracy_and_performance",
                                     "simulation_results",
                                     "primary_1.75_foursevenths_or_tablesv1.rds"))



acc_one <- read_rds(here::here("R",
                               "simulation_accuracy_and_performance",
                               "simulation_results",
                               "primary_1.75_acc_tablesv2_0.4sens.rds"))%>%
  mutate(id = 5)
summary_one <- read_rds(here::here("R",
                                   "simulation_accuracy_and_performance",
                                   "simulation_results",
                                   "primary_1.75_or_tables_0.4sens.rds"))%>%
  mutate(id = 5)


acc_three <- read_rds(here::here("R",
                                 "simulation_accuracy_and_performance",
                                 "simulation_results",
                               "primary_1.75_acc_tablesv2_0.5sens.rds"))%>%
  mutate(id = 5)
summary_three <- read_rds(here::here("R",
                                     "simulation_accuracy_and_performance",
                                     "simulation_results",
                                   "primary_1.75_or_tables_0.5sens.rds"))%>%
  mutate(id = 5)

acc_four <- read_rds(here::here("R",
                                "simulation_accuracy_and_performance",
                                "simulation_results",
                               "primary_1.75_acc_tablesv2_0.7sens.rds"))%>%
  mutate(id = 5)
summary_four <- read_rds(here::here("R",
                                    "simulation_accuracy_and_performance",
                                    "simulation_results", 
                                   "primary_1.75_or_tables_0.7sens.rds"))%>%
  mutate(id = 5)



summary <-summary_twov1%>%           
  bind_rows(.id = "id") %>%
  rbind(summary_one, summary_three, summary_four) %>%
  filter(setting != "primary_foursevenths")

acc <- acc_twov1 %>%
  bind_rows(.id = "id") %>%
  rbind(acc_one, acc_three, acc_four) %>%
  filter(setting != "primary_foursevenths")

summary <- summary %>%
  arrange(setting) %>%
  mutate(percent_bias = mean_bias/ true_val) %>%
  dplyr::select(-mean_lb, 
                - mean_ub,
                -mean_diff, 
                - id, 
                - percent_reject_low, 
                - percent_reject_high,
                -mean_bias) %>%
  dplyr::select(mean_OR, 
                percent_bias, 
                sq_MSE, 
                percent_reject,
                mean_ci_width,
                percent_coverage,
                mean_num_samples,
                num_trials,
                type, 
                setting,
                true_val
  ) %>%
  filter(type == "transphylo")

xtable(acc)
xtable(summary)


# acc graph ---------------------------------------------------------------


acc_graph_data <- acc 

acc_graph_data$true_infector <- as.character(acc_graph_data$true_infector)
acc_graph_data$true_infector[acc_graph_data$true_infector == "TRUE"] <- "Sensitivity"
acc_graph_data$true_infector[acc_graph_data$true_infector == "FALSE"] <- "Specificity"

acc_graph_data <- acc_graph_data %>%
  rename("Metric" = "true_infector")



acc_graph_data$setting[acc_graph_data$setting == "primary_1.75"] <- "Default (0.6)"
acc_graph_data$setting[acc_graph_data$setting == "16active_1.75diff_0.4"] <- "0.4"
acc_graph_data$setting[acc_graph_data$setting == "16active_1.75diff_0.5"] <- "0.5"
acc_graph_data$setting[acc_graph_data$setting == "16active_1.75diff_0.7"] <- "0.7"

level_list <- c("0.4", "0.5", "Default (0.6)", "0.7")

acc_graph_data$setting <- as.factor(acc_graph_data$setting)
acc_graph_data$setting <- factor(acc_graph_data$setting, levels = level_list)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


acc_graph <- acc_graph_data %>%
  rename("Cutoff" = "setting") %>%
  ggplot(aes(x = Metric, y = mean_percent, color = Cutoff)) + 
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  theme_bw() +
  ylim(c(0,1)) +
  ggtitle("Sensitivity and specificity with varying label cutoffs") +
  ylab("") +
  theme(text = element_text(size=15)) +
  scale_colour_manual(values = cbbPalette)


acc_graph

summary <- summary %>% dplyr::select(percent_coverage, percent_reject, percent_bias, mean_ci_width, type, setting)

