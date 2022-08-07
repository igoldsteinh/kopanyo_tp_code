# Produces visualizations of simulation results
# Based on summaries generated in summarise_simulation_results.R
library(tidyverse)
library(xtable)
library(patchwork)


# read in primary simulation results --------------------------------------
# from 16active_8year_samediff_finalcalc.R
acc_one <- read_rds(here::here("R",
                               "simulation_accuracy_and_performance",
                               "simulation_results",
                               "16active_8year_samediff_acc_tables.rds"))
summary_one <- read_rds(here::here("R",
                                   "simulation_accuracy_and_performance",
                                   "simulation_results",
                                   "16active_8year_samediff_sim_or_tables.rds"))

# from 16active_8year_1.75onlysim_finalcalc.R
acc_1.75 <- read_rds(here::here("R",
                                "simulation_accuracy_and_performance",
                                "simulation_results",
                                "16active_8year_1.75only_acc_tables.rds"))
summary_1.75 <- read_rds(here::here("R",
                                    "simulation_accuracy_and_performance",
                                    "simulation_results",
                                    "16active_8year_1.75only_sim_or_tables.rds"))

# from 16active_8year_0.3sim_finalcalc.R
acc_0.3 <- read_rds(here::here("R",
                               "simulation_accuracy_and_performance",
                               "simulation_results", 
                               "16active_8year_0.3diff_acc_tables.rds"))
summary_0.3 <- read_rds(here::here("R",
                                   "simulation_accuracy_and_performance",
                                   "simulation_results",
                                   "16active_8year_0.3diff_sim_or_tables.rds"))

# from 16active_8year_0.6sim_finalcalc.R
acc_0.6 <- read_rds(here::here("R",
                               "simulation_accuracy_and_performance",
                               "simulation_results", 
                               "16active_8year_0.6diff_acc_tables.rds"))
summary_0.6 <- read_rds(here::here("R",
                                   "simulation_accuracy_and_performance",
                                   "simulation_results",
                                   "16active_8year_0.6diff_sim_or_tables.rds"))

# from 16active_8year_3sim_finalcalc.R
acc_three <- read_rds(here::here("R",
                                 "simulation_accuracy_and_performance",
                                 "simulation_results", 
                                 "16active_8year_3diff_acc_tablesv2.rds"))
summary_three <- read_rds(here::here("R",
                                     "simulation_accuracy_and_performance",
                                     "simulation_results",
                                     "16active_8year_3diff_sim_or_tablesv2.rds")) 

summary <- summary_one %>%
  rbind(summary_three) %>%
  rbind(summary_1.75) %>%
  rbind(summary_0.3) %>%
  rbind(summary_0.6)

acc <- acc_one %>%
  rbind(acc_three) %>%
  rbind(acc_1.75) %>%
  rbind(acc_0.3) %>%
  rbind(acc_0.6)

summary <- summary %>%
  arrange(setting) %>%
  mutate(percent_bias = mean_bias/ true_val) %>%
  dplyr::select(-mean_lb, 
                - mean_ub,
                -mean_diff, 
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
  )

truth_summary <- summary %>% filter(type == "truth")
tp_summary <- summary %>% filter(type == "transphylo")
samba_summary <- summary %>% filter(type == "SAMBA Calc")
stan_summary <- summary %>% filter(type == "Stan")
# read in secondary simulation results ------------------------------------

# these come from double_density_finalcalc.R
acc_doubledens <- read_rds(here::here("R",
                                      "simulation_accuracy_and_performance",
                                      "simulation_results",
                                      "double_density_acc_tables.rds"))


summary_doubledens <- read_rds(here::here("R",
                                          "simulation_accuracy_and_performance",
                                          "simulation_results",
                                          "double_density_sim_or_tables.rds"))

# dubclusters comes from dubclusters_finalcalc.R
acc_dubclusters <- read_rds(here::here("R",
                                       "simulation_accuracy_and_performance",
                                       "simulation_results",
                                       "dubclusters_acc_tables.rds"))


summary_dubclusters <- read_rds(here::here("R",
                                           "simulation_accuracy_and_performance",
                                           "simulation_results",
                                           "dubclusers_sim_or_tables.rds"))

# 7yrsamp_1.75sim_acc_tablev1.R
acc_7yr <- read_rds(here::here("R",
                               "simulation_accuracy_and_performance",
                               "simulation_results",
                               "7yrsamp_1.75diff_acc_tablesv1.rds"))


summary_7yr <- read_rds(here::here("R",
                                   "simulation_accuracy_and_performance",
                                   "simulation_results",
                                   "7yrsamp_1.75diff_sim_or_tablesv1.rds")) 

# quad density from quadruple_density_finalcalc.R
acc_quaddens <- read_rds(here::here("R",
                                    "simulation_accuracy_and_performance",
                                    "simulation_results",
                                    "quadruple_density_acc_tables.rds"))


summary_quaddens <- read_rds(here::here("R",
                                        "simulation_accuracy_and_performance",
                                        "simulation_results", 
                                        "quadruple_density_or_tables.rds"))

summary_7yr$setting <- as.character(summary_7yr$setting)
summary_7yr$num_trials <- as.numeric(summary_7yr$num_trials)
altsummary <- data.frame(summary_quaddens[[1]]) %>%
  rbind(summary_dubclusters) %>%
  rbind(summary_doubledens) %>%
  rbind(summary_7yr)

altsummary$setting <- as.character(altsummary$setting)

altacc <- c(acc_7yr, acc_quaddens) %>%
  bind_rows(.id = "id") %>%
  dplyr::select(-id) %>%
  rbind(acc_dubclusters, acc_doubledens)

altsummary <- altsummary %>%
  arrange(setting) %>%           
  mutate(percent_bias = mean_bias/ true_val) %>%
  dplyr::select(-mean_lb, 
                - mean_ub,
                -mean_diff, 
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
  ) 

# Figure: "Operating characteristics of statistical pipelines in primary simulation settings"
# need 3 columns, value, metric, method
graph_summary <- summary %>%
  dplyr::select(type, 
                setting, 
                percent_bias, 
                percent_reject, 
                mean_ci_width, 
                percent_coverage) %>%
  pivot_longer(all_of(c("percent_bias", 
                        "percent_reject",
                        "mean_ci_width",
                        "percent_coverage")), names_to = "metric")

graph_summary$metric <- as.factor(graph_summary$metric)
level_list <- c("percent_coverage", "percent_reject", "percent_bias", "mean_ci_width")

label_list <- c("Coverage", "Reject", "Percent.Bias", "MCIW")
graph_summary$metric <- factor(graph_summary$metric, levels=level_list, labels=label_list)


graph_summary$setting[graph_summary$setting == "16active_8year_1.75diff"] <- "1.75"
graph_summary$setting[graph_summary$setting == "16active_8year_3diff"] <- "3"
graph_summary$setting[graph_summary$setting == "16active_8year_samediff"] <- "1"
graph_summary$setting[graph_summary$setting == "16active_8year_0.6diff"] <- "0.57"
graph_summary$setting[graph_summary$setting == "16active_8year_0.3diff"] <- "0.3"

level_list <- c("0.3", "0.57", "1", "1.75", "3")

graph_summary$setting <- as.factor(graph_summary$setting)
graph_summary$setting <- factor(graph_summary$setting, levels = level_list)

graph_summary$type[graph_summary$type == "SAMBA Calc"] <- "TP + ME1"
graph_summary$type[graph_summary$type == "truth"] <- "Truth + GLM"
graph_summary$type[graph_summary$type == "transphylo"] <- "TP + GLM"
graph_summary$type[graph_summary$type == "Stan"] <- "TP + ME2"

graph_summary$type <- as.factor(graph_summary$type)


level_list <- c("Truth + GLM", "TP + GLM", "TP + ME1", "TP + ME2")


graph_summary$type <- factor(graph_summary$type, levels = level_list)



cov_plot <- graph_summary %>%
            filter(metric == "Coverage") %>%
            ggplot(aes(x = type, y = value, shape = setting)) +
            geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
            geom_hline(yintercept=0.95, linetype="dashed") +
            theme_bw() +
            theme(legend.position = "none") +
            xlab("") +
            ylab("Coverage") +
            ggtitle("Coverage") +
            ylim(c(0,1)) +
            scale_shape_manual(values = c(0, 1, 2, 3, 4)) 
cov_plot


reject_plot <- graph_summary %>% 
  rename("Setting" = "setting") %>%
  filter(metric == "Reject") %>%
  ggplot(aes(x = type, y = value, shape = Setting)) +
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  geom_hline(yintercept=0.95, linetype="dashed") +
  geom_hline(yintercept=0.05, linetype="dashed") +
  theme_bw() +
  xlab("") +
  ylab("Prop. Reject") +
  ggtitle("Prop. Reject") +
scale_shape_manual(values = c(0, 1, 2, 3, 4)) 

reject_plot



bias_plot <- graph_summary %>%
  filter(metric == "Percent.Bias") %>%
  mutate(value = value * 100) %>%
  ggplot(aes(x = type, y = value, shape = setting)) +
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Method") +
  ylab("Percent Bias") +
  ggtitle("Percent Bias") +
  scale_shape_manual(values = c(0, 1, 2, 3, 4)) 

bias_plot

MCIW_plot <- graph_summary %>%
  filter(metric == "MCIW") %>%
  ggplot(aes(x = type, y = value, shape = setting)) +
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Method") +
  ylab("MCIW") +
  ggtitle("MCIW") +
  scale_shape_manual(values = c(0, 1, 2, 3, 4)) 

MCIW_plot


primary_metric_plot <- cov_plot + reject_plot + bias_plot + MCIW_plot
primary_metric_plot


ggsave(here::here("R", "simulation_accuracy_and_performance", "primary_metric_plot.pdf"), 
       plot = primary_metric_plot, 
       width = 8, 
       height = 7 )

# Figure: "Operating characteristics of statistical pipelines in secondary simulation settings"


altgraph_summary <- altsummary %>%
  dplyr::select(type, 
                setting, 
                percent_bias, 
                percent_reject, 
                mean_ci_width, 
                percent_coverage) %>%
  pivot_longer(all_of(c("percent_bias", 
                        "percent_reject",
                        "mean_ci_width",
                        "percent_coverage")), names_to = "metric")

altgraph_summary$metric <- as.factor(altgraph_summary$metric)
level_list <- c("percent_coverage", "percent_reject", "percent_bias", "mean_ci_width")

label_list <- c("Coverage", "Reject", "Percent.Bias", "MCIW")
altgraph_summary$metric <- factor(altgraph_summary$metric, levels=level_list, labels=label_list)



altgraph_summary$setting[altgraph_summary$setting == "double_density"] <- "Double Sampling Density"
altgraph_summary$setting[altgraph_summary$setting == "7yr samp"] <- "Increase Sampling Window"
altgraph_summary$setting[altgraph_summary$setting == "dubclusters"] <- "Increase Sample Size"
altgraph_summary$setting[altgraph_summary$setting == "quadruple_density"] <- "Quad Sampling Density"


altgraph_summary$type[altgraph_summary$type == "SAMBA Calc"] <- "TP + ME1"
altgraph_summary$type[altgraph_summary$type == "truth"] <- "Truth + GLM"
altgraph_summary$type[altgraph_summary$type == "transphylo"] <- "TP + GLM"
altgraph_summary$type[altgraph_summary$type == "Stan"] <- "TP + ME2"

altgraph_summary$type <- as.factor(altgraph_summary$type)


level_list <- c("Truth + GLM", "TP + GLM", "TP + ME1", "TP + ME2")


altgraph_summary$type <- factor(altgraph_summary$type, levels = level_list)

comparison <- graph_summary %>%
              filter(setting == "1.75")

comparison$setting <- as.character(comparison$setting)

altgraph_summary <- rbind(comparison, altgraph_summary)

altgraph_summary$setting[altgraph_summary$setting == "1.75"] <- "Default"

level_list <- c("Default",  
                "Double Sampling Density", 
                "Quad Sampling Density",
                "Increase Sample Size", 
                "Increase Sampling Window")

altgraph_summary$setting <- as.factor(altgraph_summary$setting)
altgraph_summary$setting <- factor(altgraph_summary$setting, levels = level_list)

# redoing this with separate bar plots so we can add in dotted lin --------

altcov_plot <- altgraph_summary %>%
  filter(metric == "Coverage") %>%
  ggplot(aes(x = type, y = value, shape = setting)) +
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  geom_hline(yintercept=0.95, linetype="dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Coverage") +
  ggtitle("Coverage") +
  ylim(c(0,1)) +
  scale_shape_manual(values = c(3, 1, 2, 4, 5)) 
altcov_plot


altreject_plot <- altgraph_summary %>% 
  rename("Setting" = "setting") %>%
  filter(metric == "Reject") %>%
  ggplot(aes(x = type, y = value, shape = Setting)) +
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  geom_hline(yintercept=0.95, linetype="dashed") +
  geom_hline(yintercept=0.05, linetype="dashed") +
  theme_bw() +
  xlab("") +
  ylab("Prop. Reject") +
  ggtitle("Prop. Reject") +
  scale_shape_manual(values = c(3, 1, 2, 4, 5)) 

altreject_plot



altbias_plot <- altgraph_summary %>%
  filter(metric == "Percent.Bias") %>%
  mutate(value = value * 100) %>%
  ggplot(aes(x = type, y = value, shape = setting)) +
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Method") +
  ylab("Percent Bias") +
  ggtitle("Percent Bias") +
  scale_shape_manual(values = c(3, 1, 2, 4, 5)) 

altbias_plot

altMCIW_plot <- altgraph_summary %>%
  filter(metric == "MCIW") %>%
  ggplot(aes(x = type, y = value, shape = setting)) +
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Method") +
  ylab("MCIW") +
  ggtitle("MCIW") +
  scale_shape_manual(values = c(3, 1, 2, 4, 5)) 

altMCIW_plot


altprimary_metric_plot <- altcov_plot + altreject_plot + altbias_plot + altMCIW_plot
altprimary_metric_plot


ggsave(here::here("R", "simulation_accuracy_and_performance", "alt_metric_plot.pdf"), 
       plot = altprimary_metric_plot, 
       width = 10, 
       height = 7 )

# Figure: "Point estimates of sensitivity and specificity of identifying infection sources"

# making a graph for the alt acc table ------------------------------------
altacc_graph_data <- altacc 

altacc_graph_data$true_infector <- as.character(altacc_graph_data$true_infector)
altacc_graph_data$true_infector[altacc_graph_data$true_infector == "TRUE"] <- "Sensitivity"
altacc_graph_data$true_infector[altacc_graph_data$true_infector == "FALSE"] <- "Specificity"

altacc_graph_data <- altacc_graph_data %>%
                     rename("Metric" = "true_infector")


altacc_graph_data$setting[altacc_graph_data$setting == "double_density"] <- "Double Sampling Density"
altacc_graph_data$setting[altacc_graph_data$setting == "7yr samp"] <- "Increase Sampling Window"
altacc_graph_data$setting[altacc_graph_data$setting == "dubclusters"] <- "Increase Clusters"
altacc_graph_data$setting[altacc_graph_data$setting == "quadruple_density"] <- "Quad Sampling Density"

altacc_graph_data <- altacc_graph_data %>%
  filter(setting %in% c("Double Sampling Density",
                        "Increase Sampling Window", 
                        "Increase Clusters", 
                        "Quad Sampling Density"))
# color palette from: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbbPalette <- c("#000000", 
                "#E69F00", 
                "#56B4E9", 
                "#009E73", 
                "#D55E00", 
                "#CC79A7",
                "#F0E442", 
                "#0072B2", 
                "#999999")
altacc_graph <- altacc_graph_data %>%
                rename("Setting" = "setting") %>%
                ggplot(aes(x = Metric, y = mean_percent, color = Setting)) + 
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  theme_bw() +
  ylim(c(0,1)) +
  ggtitle("Sens. and spec. (Secondary)") +
  ylab("")+
  theme(text = element_text(size=15)) +
  scale_colour_manual(values = cbbPalette[6:9])


altacc_graph

ggsave(here::here("R", "simulation_accuracy_and_performance",  "secondary_acc.pdf"), 
       plot = altacc_graph, 
       width = 8, 
       height = 4 )



# making a graph for primary acc table ------------------------------------
acc_graph_data <- acc 

acc_graph_data$true_infector <- as.character(acc_graph_data$true_infector)
acc_graph_data$true_infector[acc_graph_data$true_infector == "TRUE"] <- "Sensitivity"
acc_graph_data$true_infector[acc_graph_data$true_infector == "FALSE"] <- "Specificity"

acc_graph_data <- acc_graph_data %>%
  rename("Metric" = "true_infector")



acc_graph_data$setting[acc_graph_data$setting == "16active_8year_1.75diff"] <- "1.75"
acc_graph_data$setting[acc_graph_data$setting == "16active_8year_3diff"] <- "3"
acc_graph_data$setting[acc_graph_data$setting == "16active_8year_samediff"] <- "1"
acc_graph_data$setting[acc_graph_data$setting == "16active_8year_0.6diff"] <- "0.57"
acc_graph_data$setting[acc_graph_data$setting == "16active_8year_0.3diff"] <- "0.3"

level_list <- c("0.3", "0.57", "1", "1.75", "3")

acc_graph_data$setting <- as.factor(acc_graph_data$setting)
acc_graph_data$setting <- factor(acc_graph_data$setting, levels = level_list)



acc_graph <- acc_graph_data %>%
  rename("Setting" = "setting") %>%
  ggplot(aes(x = Metric, y = mean_percent, color = Setting)) + 
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  theme_bw() +
  ylim(c(0,1)) +
  ggtitle("Sens. and spec. (Primary)") +
  ylab("") +
  theme(text = element_text(size=15)) +
  scale_colour_manual(values = cbbPalette[1:5])


acc_graph

ggsave(here::here("R", "simulation_accuracy_and_performance",  "primary_acc.pdf"), 
       plot = acc_graph, 
       width = 8, 
       height = 4 )



