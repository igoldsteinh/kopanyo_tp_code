# This file contains code to analyze real data from a TB study of Botswana
# BEAST Trees at two different SNP cutoffs are used as input into TransPhylo
# The output of TransPhylo is used to generate infection source labels
# Odds ratios are obtained using three different methods of analysis
# Skip to line 377 for the final data analysis
library(ape)
library(tidyverse)
library(TransPhylo)
library(coda)
library(rstan)
library(tidybayes)
library(lubridate)
library(SAMBA)
library(patchwork)
library(ggtree)
set.seed(1)

# dictionaries and varfiles formed in creating_BEAST_clusters_snplevel5-12.R
# BEAST2 .xml files used to create maximum clade credibilit trees
# are in the BEAST2 folder


# SNP5 infection source label creation
# read in mcc trees -------------------------------------------------------

files <- c("mcctree_snp5_lineage1_cluster6",
           "mcctree_snp5_lineage1_cluster7",
           "mcctree_snp5_lineage2_cluster3",
           "mcctree_snp5_lineage4_cluster105",
           "mcctree_snp5_lineage4_cluster125",
           "mcctree_snp5_lineage4_cluster17",
           "mcctree_snp5_lineage4_cluster175",
           "mcctree_snp5_lineage4_cluster188",
           "mcctree_snp5_lineage4_cluster202",
           "mcctree_snp5_lineage4_cluster204",
           "mcctree_snp5_lineage4_cluster21",
           "mcctree_snp5_lineage4_cluster254",
           "mcctree_snp5_lineage4_cluster282",
           "mcctree_snp5_lineage4_cluster312",
           "mcctree_snp5_lineage4_cluster36",
           "mcctree_snp5_lineage4_cluster385",
           "mcctree_snp5_lineage4_cluster40",
           "mcctree_snp5_lineage4_cluster42",
           "mcctree_snp5_lineage4_cluster46",
           "mcctree_snp5_lineage4_cluster51",
           "mcctree_snp5_lineage4_cluster92")

time_tree_list <- map(files, ~read.nexus(here::here( "BEAST2", .x)))

dict5 <- read_csv(here::here("data", "kopanyo_cluster_dictionary_filtered.csv")) %>%
  filter(!is.na(mean_posterior_clock))

cluster_key <- read_csv(here::here("data", "snp5_sequence_cluster_crosswalk.csv"))


max_dates <- cluster_key %>% 
  group_by(cluster_name) %>%
  summarise(max_date = max(collectdt))

mcc_dates <- dict5 %>%
  left_join(max_dates, by = "cluster_name") %>%
  dplyr::select(mcc_tree_name, max_date) 



date_last_sample <- lubridate::decimal_date(mcc_dates$max_date)


p_tree_list <- map2(time_tree_list, date_last_sample,
                    ptreeFromPhylo)

# Run TransPhylo sharing both off.r and neg ---------------------------
# generation time gamma distribution parameters
w.shape = 10
w.scale = 1/10

# sampling time gamma distribution parameters
ws.shape = 10
ws.scale = 1/10

# Change mcmciterations to 100000 for paper analysis
# infer the transmission tree
res5 <- infer_multittree_share_param(p_tree_list,
                                    mcmcIterations = 100,
                                    thinning = 10,
                                    w.shape = w.shape,
                                    w.scale = w.scale,
                                    ws.shape = ws.shape,
                                    ws.scale = ws.scale,
                                    startOff.p = .5,
                                    updateOff.p = FALSE,
                                    prior_pi_a = 1,
                                    prior_pi_b = 19,
                                    share = c("neg", "off.r"),
                                    startNeg = 1.48,
                                    dateT = max(date_last_sample) + .01)

# write_rds(res, "snp5_kopanyo_archive_tp_allshare_res.rds")

# convergence diagnostics for all share -----------------------------------


get_param_estimates <- function(record, p){
  sapply(record, function(x) x[[p]])
}


mean3_full_param_results <- map(res5, ~data.frame(pi = get_param_estimates(., "pi"),
                                                 off.r = get_param_estimates(., "off.r"),
                                                 neg = get_param_estimates(., "neg"),
                                                 pTTree = get_param_estimates(., "pTTree"),
                                                 pPTree = get_param_estimates(., "pPTree"))) %>%
  map(~.x %>% mutate(log_posterior = pTTree + pPTree))

# ptree_plots <- map(mean3_full_param_results, ~plot(.$pPTree, type = "l"))
# ttree_plots <- map(mean3_full_param_results, ~plot(.$pTTree, type = "l"))
off.r_plots <- map(mean3_full_param_results, ~plot(.$off.r, type = "l"))
neg_plots <- map(mean3_full_param_results, ~plot(.$neg, type = "l"))
# pi_plots <- map(mean3_full_param_results, ~plot(.$pi, type = "l"))

lp_plots <-map(mean3_full_param_results, ~plot(.$log_posterior, type = "l"))


#pull effs
mean3_coda_res <- map(res5, convertToCoda)

mean3_eff <- map(mean3_coda_res, effectiveSize) %>%
  bind_rows(.id = "column_label")


# dat for analysis ----------------------------------------------

max_dates <- cluster_key %>% 
  group_by(cluster_name) %>%
  summarise(max_date = max(collectdt))

mcc_dates <- dict5 %>%
  left_join(max_dates, by = "cluster_name") %>%
  dplyr::select(mcc_tree_name, max_date) 



date_last_sample <- lubridate::decimal_date(mcc_dates$max_date)




names5 <- map(res5, pluck, 1, "ctree", "nam")
offspring5 <- map2(res5, names5, getOffspringDist, burnin = .5) 

prob_source5 <- map(offspring5, ~1-rowSums(. == 0)/dim(.)[2])

prob_inf_frame5 <- map2(names5, prob_source5, data.frame)%>%
  map(~.x %>% rename("names" =".x..i..", "prob_source" = ".y..i.."))%>%
  map(~.x %>% mutate(new_names = str_replace(.$names, ".*_", "_"))) %>%
  map(~.x %>% mutate(final_names = str_replace_all(.$new_names, "_", ""),
                     tp_source = prob_source > .6)) %>%
  bind_rows(.id = "id")



meta_name <- "Botswana_1426_good_quality_resistance_metadata.csv"

meta_data <- read_csv(here::here("data", meta_name)) %>%
  filter(gMixture < 2) %>%
  filter(!is.na(collectdt))


hiv_status <- meta_data %>%
  dplyr::select(SampleID, hivfinal_new)

id_crosswalk <- cluster_key %>%
  dplyr::select(seq.name, SampleID)

prob_inf_frame5 <- prob_inf_frame5 %>%
  left_join(id_crosswalk, by= c("names" = "seq.name")) %>%
  left_join(hiv_status, by = "SampleID")

#add a new hiv variable
# zero means hiv positive

dat_for_analysis5 <- prob_inf_frame5 %>%
  mutate(my_hiv = hivfinal_new - 1) %>%
  filter(!is.na(my_hiv))

dat_for_vis5 <- prob_inf_frame5 %>%
  mutate(my_hiv = hivfinal_new - 1) 

# write_csv(dat_for_analysis5, here::here("code", "R Code", "dat_for_analysis_snp5_allshare.csv"))

# write_csv(dat_for_vis5, here::here("code", "R Code", "dat_forvis_snp5_allshare.csv"))

# SNP10 infection source label creation
# read in mcctrees --------------------------------------------------------
dict10 <-read_csv(here::here("data", "snp_10_kopanyo_cluster_dictionary.csv"))

mcc_name <- str_c("mcctree", dict10$file_name, sep = "_")

mcc_file_list <- str_replace(mcc_name,"\\..*","") 
time_tree_list <- map(mcc_file_list, ~read.nexus(here::here("BEAST2", .x)))

# extract sample dates by group
# fine max date in each sample
sample_names <- map(time_tree_list, pluck, "tip.label") %>%
                map(data.frame) %>%
                bind_rows(.id = "tree") %>%
                rename("sample_name" = ".x..i..")


file_name <- data.frame(tree = 1:46, file_name = dict10$file_name)

sample_names$tree <- as.numeric(sample_names$tree)

sample_names <- sample_names %>%
                left_join(file_name, by = "tree") %>%
                mutate(date = str_sub(sample_name, start = -10))

sample_names$date <- as.Date(sample_names$date)

max_dates <- sample_names %>%
             group_by(file_name) %>%
             summarise(max_date = max(date))

date_last_sample <- lubridate::decimal_date(max_dates$max_date)

p_tree_list <- map2(time_tree_list, date_last_sample,
                    ptreeFromPhylo)



# analysis sharing both off.r and neg ---------------------------
# rescale these to years
# generation time gamma distribution parameters
w.shape = 10
w.scale = 1/10

# sampling time gamma distribution parameters
ws.shape = 10
ws.scale = 1/10

# infer the transmission tree
# for real analysis, change mcmcIterations to 100000
res10 <- infer_multittree_share_param(p_tree_list,
                                    mcmcIterations = 100,
                                    thinning = 10,
                                    w.shape = w.shape,
                                    w.scale = w.scale,
                                    ws.shape = ws.shape,
                                    ws.scale = ws.scale,
                                    startOff.p = .5,
                                    updateOff.p = FALSE,
                                    prior_pi_a = 1,
                                    prior_pi_b = 19,
                                    share = c("neg", "off.r"),
                                    startNeg = 1.48,
                                    dateT = max(date_last_sample) + .01)

# write_rds(res, "snp10_kopanyo_archive_tp_allshare_res.rds")


# convergence diagnostics -----------------------------------


get_param_estimates <- function(record, p){
  sapply(record, function(x) x[[p]])
}


mean3_full_param_results10 <- map(res10, ~data.frame(pi = get_param_estimates(., "pi"),
                                                 off.r = get_param_estimates(., "off.r"),
                                                 neg = get_param_estimates(., "neg"),
                                                 pTTree = get_param_estimates(., "pTTree"),
                                                 pPTree = get_param_estimates(., "pPTree"))) %>%
  map(~.x %>% mutate(log_posterior = pTTree + pPTree))

# ptree_plots <- map(mean3_full_param_results, ~plot(.$pPTree, type = "l"))
# ttree_plots <- map(mean3_full_param_results, ~plot(.$pTTree, type = "l"))
off.r_plots <- map(mean3_full_param_results10, ~plot(.$off.r, type = "l"))
neg_plots <- map(mean3_full_param_results10, ~plot(.$neg, type = "l"))
# pi_plots <- map(mean3_full_param_results10, ~plot(.$pi, type = "l"))

lp_plots <-map(mean3_full_param_results10, ~plot(.$log_posterior, type = "l"))


#pull effs
mean3_coda_res <- map(res10, convertToCoda)

mean3_eff <- map(mean3_coda_res, effectiveSize) %>%
  bind_rows(.id = "column_label")


# create dat for analysis ---------------------------------------


dict10 <-read_csv(here::here("data", "snp_10_kopanyo_cluster_dictionary.csv"))

mcc_name <- str_c("mcctree", dict10$file_name, sep = "_")

mcc_file_list <- str_replace(mcc_name,"\\..*","") 
time_tree_list <- map(mcc_file_list, ~read.nexus(here::here("code", "BEAST2", .x)))

# extract sample dates by group
# find max date in each sample
sample_names <- map(time_tree_list, pluck, "tip.label") %>%
  map(data.frame) %>%
  bind_rows(.id = "tree") %>%
  rename("sample_name" = ".x..i..")


file_name <- data.frame(tree = 1:46, file_name = dict10$file_name)

sample_names$tree <- as.numeric(sample_names$tree)

sample_names <- sample_names %>%
  left_join(file_name, by = "tree") %>%
  mutate(date = str_sub(sample_name, start = -10))

sample_names$date <- as.Date(sample_names$date)

max_dates <- sample_names %>%
  group_by(file_name) %>%
  summarise(max_date = max(date))

date_last_sample <- lubridate::decimal_date(max_dates$max_date)

p_tree_list <- map2(time_tree_list, date_last_sample,
                    ptreeFromPhylo)



names10 <- map(res10, pluck, 1, "ctree", "nam")
offspring10 <- map2(res10, names10, getOffspringDist, burnin = .5) 

prob_source10 <- map(offspring10, ~1-rowSums(. == 0)/dim(.)[2])

prob_inf_frame10 <- map2(names10, prob_source10, data.frame)%>%
  map(~.x %>% rename("names" =".x..i..", "prob_source" = ".y..i.."))%>%
  map(~.x %>% mutate(new_names = str_replace(.$names, ".*_", "_"))) %>%
  map(~.x %>% mutate(final_names = str_replace_all(.$new_names, "_", ""),
                     tp_source = prob_source > .6)) %>%
  bind_rows(.id = "id")

hist(prob_inf_frame10$prob_source)


meta_name <- "Botswana_1426_good_quality_resistance_metadata.csv"

meta_data <- read_csv(here::here("data", meta_name)) %>%
  filter(gMixture < 2) %>%
  filter(!is.na(collectdt))


hiv_status <- meta_data %>%
  dplyr::select(SampleID, hivfinal_new)

id_crosswalk <- sample_names %>%
  mutate(SampleID = str_sub(sample_name, end = -12)) %>%
  rename("seq.name" = "sample_name") %>%
  dplyr::select(seq.name, SampleID)

prob_inf_frame10 <- prob_inf_frame10 %>%
  left_join(id_crosswalk, by= c("names" = "seq.name")) %>%
  left_join(hiv_status, by = "SampleID")

#add a new hiv variable
# zero means hiv positive

dat_for_analysis10 <- prob_inf_frame10 %>%
  mutate(my_hiv = hivfinal_new - 1) %>%
  filter(!is.na(my_hiv))

dat_for_vis10 <- prob_inf_frame10 %>%
  mutate(my_hiv = hivfinal_new - 1) 



# write_csv(dat_for_analysis10, here::here("code", "R Code", "dat_for_analysis_snp10_allshare.csv"))
# write_csv(dat_for_vis10, here::here("code", "R Code", "dat_forvis_snp10_allshare.csv"))



# Analyze the data --------------------------------------------------------
dat_for_analysis5 <- read_csv(here::here("R",
                                         "real_data_analysis",
                                         "dat_for_analysis_snp5_allshare.csv"))

dat_for_analysis10 <- read_csv(here::here("R",
                              "real_data_analysis",
                              "dat_for_analysis_snp10_allshare.csv"))

# num infection sources
sum(dat_for_analysis5$tp_source)
sum(dat_for_analysis10$tp_source)

#test how number changes with different probability cutoffs
dat_for_analysis5 <- dat_for_analysis5 %>%
                     mutate(lowcutoff = prob_source > 0.4,
                            highcutoff = prob_source > 0.8)

dat_for_analysis10 <- dat_for_analysis10 %>%
                      mutate(lowcutoff = prob_source > 0.4, 
                             highcutoff = prob_source > 0.8)

sum(dat_for_analysis5$lowcutoff)
sum(dat_for_analysis5$highcutoff)
sum(dat_for_analysis10$lowcutoff)
sum(dat_for_analysis10$highcutoff)
# num clusters, mean counts, mean tree heights
summary10 <- dict10 %>%
  summarise(
    num_cluster = n(),
    mean_count = mean(count),
    mean_height = mean(mean_tree_height)
  )

summary5 <- dict5 %>%
  summarise(
    num_cluster = n(),
    mean_count = mean(count),
    mean_height = mean(mean_tree_height)
  )

# unadjusted tp glm 
glm5 <- glm(tp_source ~ my_hiv, data = dat_for_analysis5, family = binomial)
OR_ci5 <- exp(confint.default(glm5))
OR_ci5[2,2]-OR_ci5[2,1]

glm10 <- glm(tp_source ~ my_hiv, data = dat_for_analysis10, family = binomial)
OR_ci10 <- exp(confint.default(glm10))
OR_ci10[2,2]-OR_ci10[2,1]

# glm with covariates
# vars to include 
vars <- meta_data %>%
  dplyr::select(SampleID, smokef1, tb_everf3, alcohol_excess, agenew, genderf1) 

# new data
# mygender, 0 is male, 1 is female
# mysmoke 0 is yes smoke, 1 is no smoke
# mytb is 0 is yes had tb, 1 is no did not have tb
#myalcohol 0 is yes excess drinking, 1 is no
new_data5 <- dat_for_analysis5 %>%
  left_join(vars, by = "SampleID") %>%
  mutate(mygender = genderf1 -1,
         mysmoke = smokef1 -1,
         mytb = tb_everf3 -1,
         myalcohol = alcohol_excess -1)


glm5_cov <- glm(tp_source ~ my_hiv + mygender + mysmoke + mytb + myalcohol, data = new_data5, family = binomial)
OR_ci5_cov <- round(exp(confint.default(glm5_cov)), 2)
point_est5 <- round(exp(summary(glm5_cov)[["coefficients"]][,1]), 2)
names <- c("Intercept", "HIV (1 = No)", "Sex (1 = F)", "Smoke (1 = No)", "Past TB (1 = no)", "Alc (1= no)")

snp5_cov_coeffs <- data.frame(cbind(names, point_est5, OR_ci5_cov)) %>%
  rename("point_est" = 2, "CI2.5%" = 3, "CI97.5%" = 4)
rownames(snp5_cov_coeffs) <- NULL

new_data10 <- dat_for_analysis10 %>%
  left_join(vars, by = "SampleID") %>%
  mutate(mygender = genderf1 -1,
         mysmoke = smokef1 -1,
         mytb = tb_everf3 -1,
         myalcohol = alcohol_excess -1)

glm10_cov <- glm(tp_source ~ my_hiv + mygender + mysmoke + mytb + myalcohol, data = new_data10, family = binomial)
OR_ci10_cov <- round(exp(confint.default(glm10_cov)),2)
point_est10 <- round(exp(summary(glm10_cov)[["coefficients"]][,1]),2)

snp10_cov_coeffs <- data.frame(cbind(names, point_est10, OR_ci10_cov)) %>%
  rename("point_est" = 2, "CI2.5%" = 3, "CI97.5%" = 4)
rownames(snp10_cov_coeffs) <- NULL

# snp5
SAMBA_sens5 <- sensitivity(as.numeric(dat_for_analysis5$tp_source),
                           X = dat_for_analysis5$my_hiv,
                           0.69,
                           r = NULL,
                           weights = NULL) %>%
  pluck(1)


SAMBA_res5 <- approxdist(as.numeric(dat_for_analysis5$tp_source),
                         dat_for_analysis5$my_hiv,
                         SAMBA_sens5) %>%
  bind_rows(.id = "sim") %>%
  spread(sim, Z) %>%
  mutate(CI_Low = param - 1.96 * sqrt(variance),
         CI_High = param + 1.96 * sqrt(variance),
         OR = exp(param),
         OR_Low = exp(CI_Low),
         OR_High = exp(CI_High))

# snp10
SAMBA_sens10 <- sensitivity(as.numeric(dat_for_analysis10$tp_source),
                            X = dat_for_analysis10$my_hiv,
                            0.69,
                            r = NULL,
                            weights = NULL) %>%
  pluck(1)

SAMBA_res10 <- approxdist(as.numeric(dat_for_analysis10$tp_source),
                          dat_for_analysis10$my_hiv,
                          SAMBA_sens10) %>%
  bind_rows(.id = "sim") %>%
  spread(sim, Z) %>%
  mutate(CI_Low = param - 1.96 * sqrt(variance),
         CI_High = param + 1.96 * sqrt(variance),
         OR = exp(param),
         OR_Low = exp(CI_Low),
         OR_High = exp(CI_High))

# bayesian homebrew 
# snp5
model_objects5 <- list(N = dim(dat_for_analysis5)[1],
                       z = dat_for_analysis5$tp_source,
                       x = dat_for_analysis5$my_hiv,
                       specificity = .97,
                       sensitivity = .28)




stan_fit5 <- stan(file =here::here("R", "misclass_logistic_regression.stan"),
                  data = model_objects5,
                  seed = 45,
                  iter = 2000,
                  chains = 4)
rstan::traceplot(stan_fit5, pars = "lp__")

summary_beta15 <- stan_fit5 %>%
  spread_draws(beta1) %>%
  mean_qi() %>%
  mutate(
    OR = exp(beta1),
    OR_Low = exp(.lower),
    OR_high =exp(.upper)
  )

# snp10
model_objects10 <- list(N = dim(dat_for_analysis10)[1],
                        z = dat_for_analysis10$tp_source,
                        x = dat_for_analysis10$my_hiv,
                        specificity = .97,
                        sensitivity = .28)




stan_fit10 <- stan(file =here::here("R", "misclass_logistic_regression.stan"),
                   data = model_objects10,
                   seed = 45,
                   iter = 2000,
                   chains = 4)
rstan::traceplot(stan_fit10, pars = "lp__")

summary_beta110 <- stan_fit10 %>%
  spread_draws(beta1) %>%
  mean_qi() %>%
  mutate(
    OR = exp(beta1),
    OR_Low = exp(.lower),
    OR_high =exp(.upper)
  )


# Population summary table ------------------------------------------------

meta_name <- "Botswana_1426_good_quality_resistance_metadata.csv"

meta_data <- read_csv2(here::here("data", meta_name)) %>%
  filter(gMixture < 2) %>%
  filter(!is.na(collectdt)) %>%
  mutate(my_hiv = hivfinal_new - 1) %>%
  filter(!is.na(my_hiv))

table_data <- meta_data %>%
  dplyr::select(SampleID, Lineage, genderf1, agenew, my_hiv)

lineage_count <- table_data %>%
  group_by(Lineage) %>%
  count() %>%
  ungroup() %>%
  mutate(total = sum(n),
         percent = n/total)

gender_count <- table_data %>%
  group_by(genderf1) %>%
  count() %>%
  ungroup() %>%
  mutate(total = sum(n),
         percent = n/total)


table_data$agenew <- as.numeric(table_data$agenew)

table_data$agecat <- "Under 18"
table_data$agecat[table_data$agenew > 18 & table_data$agenew <= 35] <- "18 to 35"
table_data$agecat[table_data$agenew > 35 & table_data$agenew <= 55] <- "35 to 55"
table_data$agecat[table_data$agenew > 55 & table_data$agenew <= 75] <- "55 to 75"
table_data$agecat[table_data$agenew >75] <- "75+"

age_count <- table_data %>%
  group_by(agecat) %>%
  count() %>%
  ungroup() %>%
  mutate(total = sum(n),
         percent = n/total)


hiv_count <- table_data %>%
  group_by(my_hiv) %>%
  count() %>%
  ungroup() %>%
  mutate(total = sum(n),
         percent = n/total)

# BEAST Tree visualizations -----------------------------------------------
dat_for_vis10 <- read_csv(here::here("R", 
                                    "real_data_analysis",
                                    "dat_forvis_snp10_allshare.csv"))

dict10 <-read_csv(here::here("data", "snp_10_kopanyo_cluster_dictionary.csv"))

mcc_name <- str_c("mcctree", dict10$file_name, sep = "_")

mcc_file_list <- str_replace(mcc_name,"\\..*","") 
time_tree_list <- map(mcc_file_list, ~read.nexus(here::here("BEAST2", .x)))

# create sample names by tree
sample_names <- map(time_tree_list, pluck, "tip.label") %>%
  map(data.frame) %>%
  bind_rows(.id = "tree") %>%
  rename("sample_name" = ".x..i..")


file_name <- data.frame(tree = 1:46, file_name = dict10$file_name)

sample_names$tree <- as.numeric(sample_names$tree)

sample_names <- sample_names %>%
  left_join(file_name, by = "tree") %>%
  mutate(date = str_sub(sample_name, start = -10))

sample_names$date <- as.Date(sample_names$date)


i <- 1
for (i in 1:length(time_tree_list)) {
  time_tree_list[[i]][["tip.label"]] <- gsub("_.*","",  time_tree_list[[i]][["tip.label"]])
  time_tree_list[[i]][["tip.label"]] <- gsub(".*-","",  time_tree_list[[i]][["tip.label"]])
  
}

max_dates <- sample_names %>%
  group_by(file_name) %>%
  summarise(max_date = max(date),
            max_decimal_date = lubridate::decimal_date(max_date))


# now figure out which trees have tp infection sources in them 
is <- dat_for_vis10 %>%
  dplyr::select(names, tp_source, my_hiv)
sample_names <- sample_names %>%
  left_join(is, by = c("sample_name" = "names")) %>%
  group_by(tree) %>%
  mutate(has_source = any(tp_source == TRUE)) %>%
  left_join(max_dates, by = "file_name")

sample_names$sample_name <- gsub("_.*","",sample_names$sample_name)
sample_names$sample_name <- gsub(".*-","",sample_names$sample_name)

tree_summary <- sample_names %>% 
  group_by(tree) %>%
  summarise(
    has_source = sum(has_source) > 0,
    num_samples = n()
  ) 
# lets use 43,12,9,16

tree43 <- time_tree_list[[43]]
tree9 <- time_tree_list[[9]]
tree12 <- time_tree_list[[12]]
tree16 <- time_tree_list[[16]]

sample_names <- sample_names %>%
  mutate(graph_hiv = my_hiv == 0)
data43 <- sample_names %>%
  filter(tree == 43)

data9 <- sample_names %>%
  filter(tree == 9)

data12 <- sample_names %>%
  filter(tree == 12)

data16 <- sample_names %>%
  filter(tree==16)



tree_plot43<- ggtree(tree43, mrsd = unique(data43$max_date)) + 
  #geom_tiplab(size = 2) +
  ggtitle("SNP 10") + 
  theme_tree2() + 
  xlab("Time")


data43 <- data43 %>%
  mutate(sample_name = fct_relevel(sample_name, rev(get_taxa_name(tree_plot43))))

level43 <- get_taxa_name(tree_plot43)

hiv_status43 <-ggplot(data = data43, aes( x=1, y=sample_name)) + 
  geom_tile(aes(fill=graph_hiv)) + 
  theme_tree() +
  labs(fill = "HIV Status") +
  scale_fill_manual(values = alpha(c("blue", "red"), .3)) +
  theme(legend.position = "none")





tp_source43 <- ggplot(data43, aes(x=1, y = sample_name)) +
  geom_tile(aes(fill = tp_source)) +
  theme_tree() +
  labs(fill = "Infection Source") +
  theme(legend.position = "none")




tree_plot9<- ggtree(tree9, mrsd = unique(data9$max_date)) + 
  # geom_tiplab(size = 2) +
  ggtitle("") + 
  theme_tree2() + 
  xlab("Time")


data9 <- data9 %>%
  mutate(sample_name = fct_relevel(sample_name, rev(get_taxa_name(tree_plot9))))

level9 <- get_taxa_name(tree_plot9)

hiv_status9 <-ggplot(data = data9, aes( x=1, y=sample_name)) + 
  geom_tile(aes(fill=graph_hiv)) + 
  theme_tree() +
  labs(fill = "HIV Status") +
  scale_fill_manual(values = alpha(c("blue", "red"), .3)) +
  theme(legend.position = "none")





tp_source9 <- ggplot(data9, aes(x=1, y = sample_name)) +
  geom_tile(aes(fill = tp_source)) +
  theme_tree() +
  labs(fill = "Infection Source") +
  theme(legend.position = "none")


tree_plot12<- ggtree(tree12, mrsd = unique(data12$max_date)) + 
  # geom_tiplab(size = 2) +
  ggtitle("") + 
  theme_tree2() + 
  xlab("Time")


data12 <- data12 %>%
  mutate(sample_name = fct_relevel(sample_name, rev(get_taxa_name(tree_plot12))))

level12 <- get_taxa_name(tree_plot12)

hiv_status12 <-ggplot(data = data12, aes( x=1, y=sample_name)) + 
  geom_tile(aes(fill=graph_hiv)) + 
  theme_tree() +
  labs(fill = "HIV Status") +
  scale_fill_manual(values = alpha(c("blue", "red"), .3)) +
  theme(legend.position = "none")





tp_source12 <- ggplot(data12, aes(x=1, y = sample_name)) +
  geom_tile(aes(fill = tp_source)) +
  theme_tree() +
  labs(fill = "Infection Source") +
  theme(legend.position = "none")


tree_plot16<- ggtree(tree16, mrsd = unique(data16$max_date)) + 
  geom_tiplab(size = 2.5) +
  ggtitle("") + 
  theme_tree2() + 
  xlab("Time")


data16 <- data16 %>%
  mutate(sample_name = fct_relevel(sample_name, rev(get_taxa_name(tree_plot16))))

level16 <- get_taxa_name(tree_plot16)

hiv_status16 <-ggplot(data = data16, aes( x=1, y=sample_name)) + 
  geom_tile(aes(fill=graph_hiv)) + 
  theme_tree() +
  labs(fill = "HIV Status") +
  scale_fill_manual(values = alpha(c("blue", "red"), .3)) +
  theme(legend.position = "none")





tp_source16 <- ggplot(data16, aes(x=1, y = sample_name)) +
  geom_tile(aes(fill = tp_source)) +
  theme_tree() +
  labs(fill = "Infection Source") +
  theme(legend.position = "none")


test <- (tree_plot43 + hiv_status43 + tp_source43 +plot_layout(widths = c(19,1,1)) ) /
  (tree_plot12 + hiv_status12 + tp_source12 + plot_layout(widths = c(19,1,1))) /
  (tree_plot9 + hiv_status9 + tp_source9 + plot_layout(widths = c(19,1,1))) / 
  (tree_plot16 + hiv_status16 + tp_source16 + plot_layout(widths = c(19,1,1))) 



# repeat but for 5
files <- c("mcctree_snp5_lineage1_cluster6",
           "mcctree_snp5_lineage1_cluster7",
           "mcctree_snp5_lineage2_cluster3",
           "mcctree_snp5_lineage4_cluster105",
           "mcctree_snp5_lineage4_cluster125",
           "mcctree_snp5_lineage4_cluster17",
           "mcctree_snp5_lineage4_cluster175",
           "mcctree_snp5_lineage4_cluster188",
           "mcctree_snp5_lineage4_cluster202",
           "mcctree_snp5_lineage4_cluster204",
           "mcctree_snp5_lineage4_cluster21",
           "mcctree_snp5_lineage4_cluster254",
           "mcctree_snp5_lineage4_cluster282",
           "mcctree_snp5_lineage4_cluster312",
           "mcctree_snp5_lineage4_cluster36",
           "mcctree_snp5_lineage4_cluster385",
           "mcctree_snp5_lineage4_cluster40",
           "mcctree_snp5_lineage4_cluster42",
           "mcctree_snp5_lineage4_cluster46",
           "mcctree_snp5_lineage4_cluster51",
           "mcctree_snp5_lineage4_cluster92")

time_tree_list <- map(files, ~read.nexus(here::here("BEAST2", .x)))

dat_for_vis5 <- read_csv(here::here("R", 
                                    "real_data_analysis",
                                    "dat_forvis_snp5_allshare.csv"))



# create sample names by tree
sample_names <- map(time_tree_list, pluck, "tip.label") %>%
  map(data.frame) %>%
  bind_rows(.id = "tree") %>%
  rename("sample_name" = ".x..i..")


file_name <- data.frame(tree = 1:46, file_name = dict10$file_name)

sample_names$tree <- as.numeric(sample_names$tree)

sample_names <- sample_names %>%
  left_join(file_name, by = "tree") %>%
  mutate(date = str_sub(sample_name, start = -10))

sample_names$date <- as.Date(sample_names$date)


i <- 1
for (i in 1:length(time_tree_list)) {
  time_tree_list[[i]][["tip.label"]] <- gsub("_.*","",  time_tree_list[[i]][["tip.label"]])
  time_tree_list[[i]][["tip.label"]] <- gsub(".*-","",  time_tree_list[[i]][["tip.label"]])
  
}

max_dates <- sample_names %>%
  group_by(file_name) %>%
  summarise(max_date = max(date),
            max_decimal_date = lubridate::decimal_date(max_date))


# now figure out which trees have tp infection sources in them 
is <- dat_for_vis10 %>%
  dplyr::select(names, tp_source, my_hiv)
sample_names <- sample_names %>%
  left_join(is, by = c("sample_name" = "names")) %>%
  group_by(tree) %>%
  mutate(has_source = any(tp_source == TRUE)) %>%
  left_join(max_dates, by = "file_name")

sample_names$sample_name <- gsub("_.*","",sample_names$sample_name)
sample_names$sample_name <- gsub(".*-","",sample_names$sample_name)

tree_summary <- sample_names %>% 
  group_by(tree) %>%
  summarise(
    has_source = sum(has_source) > 0,
    num_samples = n()
  ) 
# lets use 43,12,9,16

tree17 <- time_tree_list[[17]]
tree11 <- time_tree_list[[11]]
tree5 <- time_tree_list[[5]]
tree12_snp5 <- time_tree_list[[12]]

sample_names <- sample_names %>%
  mutate(graph_hiv = my_hiv == 0)
data17 <- sample_names %>%
  filter(tree == 17)

data11 <- sample_names %>%
  filter(tree == 11)

data12_snp5 <- sample_names %>%
  filter(tree == 12)

data5 <- sample_names %>%
  filter(tree==5)



tree_plot17<- ggtree(tree17, mrsd = unique(data17$max_date)) + 
  # geom_tiplab(size = 2) +
  ggtitle("SNP 5") + 
  theme_tree2() + 
  xlab("Time")


data17 <- data17 %>%
  mutate(sample_name = fct_relevel(sample_name, rev(get_taxa_name(tree_plot17))))

level17 <- get_taxa_name(tree_plot17)

hiv_status17 <-ggplot(data = data17, aes( x=1, y=sample_name)) + 
  geom_tile(aes(fill=graph_hiv)) + 
  theme_tree() +
  labs(fill = "HIV Status") +
  scale_fill_manual(values = alpha(c("blue", "red"), .3)) +
  theme(legend.position = "none")





tp_source17 <- ggplot(data17, aes(x=1, y = sample_name)) +
  geom_tile(aes(fill = tp_source)) +
  theme_tree() +
  labs(fill = "Infection Source") +
  theme(legend.position = "none")




tree_plot11<- ggtree(tree11, mrsd = unique(data11$max_date)) + 
  #geom_tiplab(size = 2) +
  ggtitle("") + 
  theme_tree2() + 
  xlab("Time")


data11 <- data11 %>%
  mutate(sample_name = fct_relevel(sample_name, rev(get_taxa_name(tree_plot11))))

level11 <- get_taxa_name(tree_plot11)

hiv_status11 <-ggplot(data = data11, aes( x=1, y=sample_name)) + 
  geom_tile(aes(fill=graph_hiv)) + 
  theme_tree() +
  labs(fill = "HIV Status") +
  scale_fill_manual(values = alpha(c("blue", "red"), .3)) +
  theme(legend.position = "none")





tp_source11 <- ggplot(data11, aes(x=1, y = sample_name)) +
  geom_tile(aes(fill = tp_source)) +
  theme_tree() +
  labs(fill = "Infection Source") +
  theme(legend.position = "none")


tree_plot12_snp5<- ggtree(tree12_snp5, mrsd = unique(data12_snp5$max_date)) + 
  geom_tiplab(size = 2.5) +
  ggtitle("") + 
  theme_tree2() + 
  xlab("Time")


data12_snp5 <- data12_snp5 %>%
  mutate(sample_name = fct_relevel(sample_name, rev(get_taxa_name(tree_plot12_snp5))))

level12_snp5 <- get_taxa_name(tree_plot12_snp5)

hiv_status12_snp5 <-ggplot(data = data12_snp5, aes( x=1, y=sample_name)) + 
  geom_tile(aes(fill=graph_hiv)) + 
  theme_tree() +
  labs(fill = "HIV Status") +
  scale_fill_manual(values = alpha(c("blue", "red"), .3)) +
  theme(legend.position = "none")





tp_source12_snp5 <- ggplot(data12_snp5, aes(x=1, y = sample_name)) +
  geom_tile(aes(fill = tp_source)) +
  theme_tree() +
  labs(fill = "Infection Source") +
  theme(legend.position = "none")


tree_plot5<- ggtree(tree5, mrsd = unique(data5$max_date)) + 
  #geom_tiplab(size = 1.5) +
  ggtitle("") + 
  theme_tree2() + 
  xlab("Time")


data5 <- data5 %>%
  mutate(sample_name = fct_relevel(sample_name, rev(get_taxa_name(tree_plot5))))

level5 <- get_taxa_name(tree_plot5)

hiv_status5 <-ggplot(data = data5, aes( x=1, y=sample_name)) + 
  geom_tile(aes(fill=graph_hiv)) + 
  theme_tree() +
  labs(fill = "HIV Status") +
  scale_fill_manual(values = alpha(c("blue", "red"), .3)) 





tp_source5 <- ggplot(data5, aes(x=1, y = sample_name)) +
  geom_tile(aes(fill = tp_source)) +
  theme_tree() +
  labs(fill = "Infection Source")

test <- (tree_plot43 + hiv_status43 + tp_source43 +plot_layout(widths = c(19,1,1)) ) /
  (tree_plot12_snp5 + hiv_status12_snp5 + tp_source12_snp5 + plot_layout(widths = c(19,1,1))) /
  (tree_plot9 + hiv_status9 + tp_source9 + plot_layout(widths = c(19,1,1))) / 
  (tree_plot16 + hiv_status16 + tp_source16 + plot_layout(widths = c(19,1,1))) 

full_tree_plot <- (tree_plot43 + hiv_status43 + tp_source43 +tree_plot17 + hiv_status17 + tp_source17 +plot_layout(widths = c(19,1,1, 19, 1, 1)) ) /
  (tree_plot12 + hiv_status12 + tp_source12 +tree_plot11 + hiv_status11 + tp_source11 + plot_layout(widths = c(19,1,1, 19, 1, 1))) / 
  (tree_plot9 + hiv_status9 + tp_source9 +tree_plot5 + hiv_status5 + tp_source5 + plot_layout(widths = c(19,1,1, 19, 1, 1))) /
  (tree_plot16 + hiv_status16 + tp_source16 + tree_plot12_snp5 + hiv_status12_snp5 + tp_source12_snp5 + plot_layout(widths = c(19,1,1, 19, 1, 1))) + plot_layout(guides = "collect")

plain_tree <- ggtree(tree12_snp5, mrsd = unique(data12_snp5$max_date), size = 4) + 
  ggtitle("") + 
  theme_tree2() + 
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 32),
        axis.ticks.x=element_blank()) +
  xlab("Time")

ggsave("plain_tree.png", plain_tree, width = 5, height = 5)
ggsave("full_tree_plotv2.pdf", full_tree_plot, width = 10, height = 10)


