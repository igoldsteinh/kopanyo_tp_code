library(ape)
library(tidyverse)
library(TransPhylo)
library(coda)
library(rstan)
library(tidybayes)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


#code for running transphylo on snp level 5 kopanyo archived data
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

time_tree_list <- map(files, ~read.nexus(here::here("code", "BEAST2", .x)))

# read.nexus(here::here("code", "BEAST2", "mcctree_snp5_lineage1_cluster6"))

# ok this is sort of a thing I forgot about 
# need to find the max sample date for each cluster's tree
dict <- read_csv(here::here("data", "kopanyo_cluster_dictionary_filtered.csv")) %>%
        filter(!is.na(mean_posterior_clock))

cluster_key <- read_csv(here::here("data", "snp5_sequence_cluster_crosswalk.csv"))


max_dates <- cluster_key %>% 
             group_by(cluster_name) %>%
             summarise(max_date = max(collectdt))

mcc_dates <- dict %>%
             left_join(max_dates, by = "cluster_name") %>%
             dplyr::select(mcc_tree_name, max_date) 
             


date_last_sample <- lubridate::decimal_date(mcc_dates$max_date)


p_tree_list <- map2(time_tree_list, date_last_sample,
                    ptreeFromPhylo)



# run transphylo ----------------------------------------------------------
# rescale these to years
# generation time gamma distribution parameters
w.shape = 10
w.scale = 1/10

# sampling time gamma distribution parameters
ws.shape = 10
ws.scale = 1/10

# infer the transmission tree
res <- infer_multittree_share_param(p_tree_list,
                                    mcmcIterations = 10000,
                                    w.shape = w.shape,
                                    w.scale = w.scale,
                                    ws.shape = ws.shape,
                                    ws.scale = ws.scale,
                                    startOff.p = .5,
                                    updateOff.p = FALSE,
                                    prior_pi_a = 1,
                                    prior_pi_b = 19,
                                    share = c("neg"),
                                    startNeg = 1.48,
                                    dateT = max(date_last_sample) + .01)


write_rds(res, "snp5_kopanyo_archive_tp_res.rds")
write_rds(res, "snp5_kopanyo_archive_tp_res_nor.rds")


# diagnostics -------------------------------------------------------------
res <- read_rds("snp5_kopanyo_archive_tp_res.rds")
get_param_estimates <- function(record, p){
  sapply(record, function(x) x[[p]])
}


mean3_full_param_results <- map(res, ~data.frame(pi = get_param_estimates(., "pi"),
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
mean3_coda_res <- map(res, convertToCoda)

mean3_eff <- map(mean3_coda_res, effectiveSize) %>%
  bind_rows(.id = "column_label")

names <- map(res, pluck, 1, "ctree", "nam")
offspring <- map2(res, names, getOffspringDist, burnin = .5) 
times <-map2(res, names, getInfectionTimeDist, burnin = .5) %>%
  map(data.frame) %>%
  bind_rows(.id = "sim")

prob_source <- map(offspring, ~1-rowSums(. == 0)/dim(.)[2])

prob_inf_frame <- map2(names, prob_source, data.frame)%>%
  map(~.x %>% rename("names" =".x..i..", "prob_source" = ".y..i.."))%>%
  map(~.x %>% mutate(new_names = str_replace(.$names, ".*_", "_"))) %>%
  map(~.x %>% mutate(final_names = str_replace_all(.$new_names, "_", ""),
                     tp_source = prob_source > .6)) %>%
  bind_rows(.id = "id")

hist(prob_inf_frame$prob_source)


meta_name <- "Botswana_1426_good_quality_resistance_metadata.csv"

meta_data <- read_csv2(here::here("data", meta_name)) %>%
  filter(gMixture < 2) %>%
  filter(!is.na(collectdt))


hiv_status <- meta_data %>%
              dplyr::select(SampleID, hivfinal_new)

id_crosswalk <- cluster_key %>%
                dplyr::select(seq.name, SampleID)

prob_inf_frame <- prob_inf_frame %>%
                  left_join(id_crosswalk, by= c("names" = "seq.name")) %>%
                  left_join(hiv_status, by = "SampleID")

#add a new hiv variable
# zero means hiv positive

dat_for_analysis <- prob_inf_frame %>%
                  mutate(my_hiv = hivfinal_new - 1) %>%
                  filter(!is.na(my_hiv))


source_by_group <- dat_for_analysis %>%
                  group_by(my_hiv) %>%
                  summarise(num_source = sum(tp_source))


glm_model <- glm(tp_source ~ my_hiv, family = binomial, data = dat_for_analysis)


OR_ci <- exp(confint(glm_model))
OR_ci

model_objects <- list(N = dim(dat_for_analysis)[1],
                      z = dat_for_analysis$tp_source,
                      x = dat_for_analysis$my_hiv,
                      specificity = .90,
                      sensitivity = .32)




stan_fit <- stan(file =here::here("code", "R Code", "misclass_logistic_regression.stan"),
                                          data = model_objects,
                                          seed = 45,
                                          iter = 2000,
                                          chains = 4)
rstan::traceplot(stan_fit, pars = "lp__")

summary_beta1 <- stan_fit %>%
                 spread_draws(beta1) %>%
                 mean_qi() %>%
                 mutate(
                   OR = exp(beta1),
                   OR_Low = exp(.lower),
                   OR_high =exp(.upper)
                 )

summary_beta1

# using SAMBA -------------------------------------------------------------
# 0.6484173
dat_for_analysis5<- read_csv(here::here("code", "R Code", "dat_for_analysis_snp5.csv"))


SAMBA_sens <- sensitivity(as.numeric(dat_for_analysis5$tp_source),
                          X = dat_for_analysis5$my_hiv,
                          0.6484173,
                          r = NULL,
                          weights = NULL) %>%
  pluck(1)


SAMBA_res <- approxdist(as.numeric(dat_for_analysis5$tp_source),
                        dat_for_analysis5$my_hiv,
                        SAMBA_sens) %>%
  bind_rows(.id = "sim") %>%
  spread(sim, Z) %>%
  mutate(CI_Low = param - 1.96 * sqrt(variance),
         CI_High = param + 1.96 * sqrt(variance),
         OR = exp(param),
         OR_Low = exp(CI_Low),
         OR_High = exp(CI_High))

# repeating analysis sharing both off.r and neg ---------------------------
# rescale these to years
# generation time gamma distribution parameters
w.shape = 10
w.scale = 1/10

# sampling time gamma distribution parameters
ws.shape = 10
ws.scale = 1/10

# infer the transmission tree
res <- infer_multittree_share_param(p_tree_list,
                                    mcmcIterations = 100000,
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

write_rds(res, "snp5_kopanyo_archive_tp_allshare_res.rds")

# convergence diagnostics for all share -----------------------------------

res <- read_rds(here::here("code", "R Code", "snp5_kopanyo_archive_tp_allshare_res.rds"))

get_param_estimates <- function(record, p){
  sapply(record, function(x) x[[p]])
}


mean3_full_param_results <- map(res, ~data.frame(pi = get_param_estimates(., "pi"),
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
mean3_coda_res <- map(res, convertToCoda)

mean3_eff <- map(mean3_coda_res, effectiveSize) %>%
  bind_rows(.id = "column_label")


# dat for analysis all sahre ----------------------------------------------

dict <- read_csv(here::here("data", "kopanyo_cluster_dictionary_filtered.csv")) %>%
  filter(!is.na(mean_posterior_clock))

cluster_key <- read_csv(here::here("data", "snp5_sequence_cluster_crosswalk.csv"))


max_dates <- cluster_key %>% 
  group_by(cluster_name) %>%
  summarise(max_date = max(collectdt))

mcc_dates <- dict %>%
  left_join(max_dates, by = "cluster_name") %>%
  dplyr::select(mcc_tree_name, max_date) 



date_last_sample <- lubridate::decimal_date(mcc_dates$max_date)



res5 <- read_rds("snp5_kopanyo_archive_tp_allshare_res.rds")

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

meta_data <- read_csv2(here::here("data", meta_name)) %>%
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

write_csv(dat_for_analysis5, here::here("code", "R Code", "dat_for_analysis_snp5_allshare.csv"))

write_csv(dat_for_vis5, here::here("code", "R Code", "dat_forvis_snp5_allshare.csv"))

