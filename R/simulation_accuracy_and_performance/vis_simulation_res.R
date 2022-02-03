library(tidyverse)
library(xtable)
library(patchwork)

# from 16active_8year_diffsim_finalcalc.R
acc_one <- read_rds(here::here("code", 
                           "R Code", 
                           "16active_8year_acc_tables.rds"))
summary_one <- read_rds(here::here("code", 
                               "R Code", 
                               "16active_8year_sim_or_tables.rds"))
# from 16active_8year_1.75sim_finalcalc.R
acc_twov1 <- read_rds(here::here("code", 
                               "R Code", 
                               "16active_8year_1.75diff_acc_tablesv1.rds"))
summary_twov1 <- read_rds(here::here("code", 
                                   "R Code", 
                                   "16active_8year_1.75diff_sim_or_tablesv1.rds"))

# # coverage remains the same mean bias differs by a little bit but no takeaways are changed
# acc_twov2 <- read_rds(here::here("code", 
#                                  "R Code", 
#                                  "16active_8year_1.75diff_acc_tablesv2.rds"))
# summary_twov2 <- read_rds(here::here("code", 
#                                      "R Code", 
#                                      "16active_8year_1.75diff_sim_or_tablesv2.rds"))

# from 16active_8year_3sim_finalcalc.R
acc_three <- read_rds(here::here("code", 
                               "R Code", 
                               "16active_8year_3diff_acc_tablesv2.rds")) %>%
              mutate(id = 5)
summary_three <- read_rds(here::here("code", 
                                   "R Code", 
                                   "16active_8year_3diff_sim_or_tablesv2.rds")) %>%
                  mutate(id = 5)


summary <- c(summary_one, summary_twov1)%>%           
            bind_rows(.id = "id") %>%
            rbind(summary_three)

acc <- c(acc_one, acc_twov1) %>%
       bind_rows(.id = "id") %>%
       rbind(acc_three)

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
                         )


# trying to visualize this table as a figure  -----------------------------

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
graph_summary$setting[graph_summary$setting == "16active_8yearsim_0.6diff"] <- "4/7"
graph_summary$setting[graph_summary$setting == "16active_8yearsim_0.3diff"] <- "0.3"

level_list <- c("0.3", "4/7", "1", "1.75", "3")

graph_summary$setting <- as.factor(graph_summary$setting)
graph_summary$setting <- factor(graph_summary$setting, levels = level_list)

graph_summary$type[graph_summary$type == "SAMBA Calc"] <- "TP + ME1"
graph_summary$type[graph_summary$type == "truth"] <- "Truth + GLM"
graph_summary$type[graph_summary$type == "transphylo"] <- "TP + GLM"
graph_summary$type[graph_summary$type == "Stan"] <- "TP + ME2"

graph_summary$type <- as.factor(graph_summary$type)


level_list <- c("Truth + GLM", "TP + GLM", "TP + ME1", "TP + ME2")


graph_summary$type <- factor(graph_summary$type, levels = level_list)
bar_plot <- graph_summary %>%
            ggplot(aes(x = type, y = value, shape = setting)) +
            geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
            facet_wrap(~metric, scales = "free_y", labeller = label_parsed) +
            theme_bw()
bar_plot


# redoing this with separate bar plots so we can add in dotted lin --------

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


ggsave(here::here("code", "R Code", "primary_metric_plot.pdf"), 
       plot = primary_metric_plot, 
       width = 8, 
       height = 7 )


# let's make the same graph for the alternate settings --------------------
# these come from altsetting_1.75sim_finalcalc.R
altacc_one <- read_rds(here::here("code", 
                               "R Code", 
                               "altsetting_1.75diff_acc_tablesv1.rds"))


altsummary_one <- read_rds(here::here("code", 
                                   "R Code", 
                                   "altsetting_1.75diff_sim_or_tablesv1.rds"))

# dubclusters comes from dubclusters_1.75sim_results_tablev1.R
altacc_two <- read_rds(here::here("code", 
                               "R Code", 
                               "dubclusters_1.75diff_acc_tablesv1.rds"))


altsummary_two <- read_rds(here::here("code", 
                                   "R Code", 
                                   "dubclusters_1.75diff_sim_or_tablesv1.rds"))

# 7yrsamp_1.75sim_acc_tablev1.R
altacc_three <- read_rds(here::here("code", 
                               "R Code", 
                               "7yrsamp_1.75diff_acc_tablesv1.rds"))


altsummary_three <- read_rds(here::here("code", 
                                   "R Code", 
                                   "7yrsamp_1.75diff_sim_or_tablesv1.rds")) %>%
                    mutate(id = "three")

altsummary <- c(altsummary_one, altsummary_two) %>%
  bind_rows(.id = "id") %>%
  rbind(altsummary_three)

altacc <- c(altacc_one, altacc_two, altacc_three) %>%
  bind_rows(.id = "id") 

altsummary <- altsummary %>%
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
  filter(setting != "32active_4yearsim_1.75diff" & setting!= "16active_4yearsim_1.75diff")

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


altgraph_summary$setting[altgraph_summary$setting == "32active_8yearsim_1.75"] <- "Increase Sampling Density"
altgraph_summary$setting[altgraph_summary$setting == "7yr samp"] <- "Increase Sampling Window"
altgraph_summary$setting[altgraph_summary$setting == "16active_8yearsim_1.75diff_dubclusters"] <- "Increase Clusters"


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

level_list <- c("Default", "Increase Sampling Density", "Increase Clusters", "Increase Sampling Window")

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
  scale_shape_manual(values = c(3, 1, 2, 4)) 
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
  scale_shape_manual(values = c(3, 1, 2, 4)) 

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
  scale_shape_manual(values = c(3, 1, 2, 4)) 

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
  scale_shape_manual(values = c(3, 1, 2, 4)) 

altMCIW_plot


altprimary_metric_plot <- altcov_plot + altreject_plot + altbias_plot + altMCIW_plot
altprimary_metric_plot


ggsave(here::here("code", "R Code", "alt_metric_plot.pdf"), 
       plot = altprimary_metric_plot, 
       width = 10, 
       height = 7 )


# making a graph for the alt acc table ------------------------------------
altacc_graph_data <- altacc 

altacc_graph_data$true_infector <- as.character(altacc_graph_data$true_infector)
altacc_graph_data$true_infector[altacc_graph_data$true_infector == "TRUE"] <- "Sensitivity"
altacc_graph_data$true_infector[altacc_graph_data$true_infector == "FALSE"] <- "Specificity"

altacc_graph_data <- altacc_graph_data %>%
                     rename("Metric" = "true_infector")


altacc_graph_data$setting[altacc_graph_data$setting == "32active_8yearsim_1.75"] <- "Increase Sampling Density"
altacc_graph_data$setting[altacc_graph_data$setting == "7yr samp"] <- "Increase Sampling Window"
altacc_graph_data$setting[is.na(altacc_graph_data$setting)] <- "Increase Clusters"

altacc_graph_data <- altacc_graph_data %>%
                     filter(setting %in% c("Increase Sampling Density",
                                           "Increase Sampling Window", 
                                           "Increase Clusters"))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
altacc_graph <- altacc_graph_data %>%
                rename("Setting" = "setting") %>%
                ggplot(aes(x = Metric, y = mean_percent, color = Setting)) + 
  geom_jitter(width = 0, height = 0.005, stroke = 0.75) +
  theme_bw() +
  ylim(c(0,1)) +
  ggtitle("Sens. and spec. (Secondary)") +
  ylab("")+
  theme(text = element_text(size=15)) +
  scale_colour_manual(values = cbbPalette[5:7])


altacc_graph


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
acc_graph_data$setting[acc_graph_data$setting == "16active_8yearsim_0.6diff"] <- "4/7"
acc_graph_data$setting[acc_graph_data$setting == "16active_8yearsim_0.3diff"] <- "0.3"

level_list <- c("0.3", "4/7", "1", "1.75", "3")

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


combined_acc <- acc_graph + altacc_graph +
  plot_layout(guides = "collect")

ggsave(here::here("code", "R Code", "combined_acc.pdf"), 
       plot = combined_acc, 
       width = 10, 
       height = 4 )

