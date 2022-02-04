# This file is for splitting fasta files into cluster groupings and renaming them for input into BEAST

library(tidyverse)
library(ape)
library(transcluster)
library(phylotools)
library(treedater)
library(lubridate)
library(here)
library(aplot)
library(ggtree)
library(xtable)
set.seed(1)


# reading in data ---------------------------------------------------------


samples1 <- read.dna(file=here("data", "Botswana_lineage_1_86_samples.fasta"), 
                     format = "fasta") %>% 
  `rownames<-`(str_extract(labels(.), "^BTB-\\d+"))

samples2 <- read.dna(file=here("data", "Botswana_lineage_2_76_samples.fasta"), 
                     format = "fasta") %>% 
  `rownames<-`(str_extract(labels(.), "^BTB-\\d+"))

samples3 <- read.dna(file=here("data", "Botswana_lineage_3_12_samples.fasta"), 
                     format = "fasta") %>% 
  `rownames<-`(str_extract(labels(.), "^BTB-\\d+"))

samples4 <- read.dna(file=here("data", "Botswana_lineage_4_1252_samples.fasta"), 
                     format = "fasta") %>% 
  `rownames<-`(str_extract(labels(.), "^BTB-\\d+"))



meta_name <- "Botswana_1426_good_quality_resistance_metadata.csv"

meta_data <- read_csv(here("data", meta_name)) %>%
  filter(gMixture < 2) %>%
  filter(!is.na(collectdt))

num_samples <- dim(meta_data)[1]

samples <- list(samples1, samples2, samples3, samples4)


samples <- map(samples, ~.[rownames(.) %in% meta_data$SampleID,])

meta_data <- meta_data %>%
  mutate(lineage_rough = str_sub(Lineage_detailed,1, 1)) %>%
  split(.$lineage_rough)


# clustering --------------------------------------------------------------
# snp matrix from 2021_03_12_1426_sample_clusters_and_snps
# actual code 
# snp_matrix <- map(samples, dist.gene, method = "pairwise")%>%
# map(as.matrix)

snp_matrix <- read_rds(here::here("code", 
                                  "R Code", 
                                  str_c(as.character(num_samples), 
                                        "_samples_snp_matrix_list.rds")))
ids <- map(samples, rownames)

lambda <- 1.5
beta <- 2.0

dates <- map(meta_data, ~lubridate::year(.x$collectdt))

my_model <- pmap(list(ids, dates, snp_matrix), createModel)
# pick which SNP cutoffs to use
my_model <- map(my_model, setSNPThresholds,c(5,10))
# pick which transmission cutoffs to use
my_model <- map(my_model, setTransThresholds, seq(4, 24, by=2))
# choose clock rate (lambda) and transmission rate (beta)
my_model <- map(my_model, setParams, lambda = lambda, beta = beta) 

my_SNP_clusters <- map(my_model, makeSNPClusters, writeFile = F)
# my_trans_clusters <- map(my_model, makeTransClusters, writeFile = F)

snp_levels <- c(5, 6,7,8,9,10,11,12)
id <- meta_data$SampleID

ids <- map(samples, rownames)


id_list_snp <- map(ids, list) %>%
  map(rep, times = 8)
snp_frame <- map2(my_SNP_clusters, id_list_snp, ~map2(.x, .y, data.frame))
snp_frame <- map(snp_frame, bind_rows, .id="id")%>%
  map(~.x %>%
        rename("cluster"= ".x..i..",
               "snp_level"="id",
               "id" = ".y..i.."))


# splitting samples up into clusters at different snp levels ---------------
# I apologize to Mine but I'm lazy and I'm just going to loop it
snp_levels<- c(5, 6,7,8,9,10,11,12)

samples1 <- read.fasta(file=here("data", "Botswana_lineage_1_86_samples.fasta"))

samples2 <- read.fasta(file=here("data", "Botswana_lineage_2_76_samples.fasta"))


samples3 <- read.fasta(file=here("data", "Botswana_lineage_3_12_samples.fasta"))

samples4 <- read.fasta(file=here("data", "Botswana_lineage_4_1252_samples.fasta"))

for (i in 1: length(snp_levels)) {
  level <- snp_levels[i]
  snp_level <- str_c("snp", level, sep="")
  
  snp_frame_level <- map(snp_frame, ~.x %>%
                     filter(snp_level == level) %>%
                     group_by(cluster) %>% 
                     mutate(count = n()) %>% 
                     filter(count >= 4))
  
  samples_data_frame <- list(samples1, samples2, samples3, samples4)
  
  
  id_times <- map(meta_data, ~.x %>%
                    dplyr::select(SampleID, collectdt, lineage_rough) %>%
                    mutate(date_id = str_c(SampleID, as.character(collectdt), sep = "_")))
  
  
  
  samples_data_frame <- pmap(list(samples_data_frame, id_times, snp_frame_level), ~..1 %>% 
                               left_join(..2, by = c("seq.name" = "SampleID")) %>%
                               filter(seq.name %in% ..3$id) %>%
                               dplyr::select(date_id, seq.text) %>%
                               rename("seq.name" = "date_id"))
  
  
  
  cluster_names <- map2(snp_frame_level, id_times, ~.x %>%
                          dplyr::select(id, cluster) %>%
                          left_join(.y, by = c("id" = "SampleID")) %>%
                          rename("seq.name" = "date_id") %>%
                          mutate(cluster_name = str_c(snp_level,
                                                      "lineage", 
                                                      lineage_rough, 
                                                      "cluster", 
                                                      cluster, 
                                                      sep = "_")) %>%
                          ungroup() %>%
                          dplyr::select(seq.name, cluster_name))
  
  map2(samples_data_frame, cluster_names, split_dat)
  
  match <- str_c("^snp", level, sep="")
  file_list <- list.files("C:/Users/fiddl/Documents/kopanyo-archived-phylo/", pattern = match)
  new_samples <- map(file_list, ~read.dna(here::here(.x), format = "fasta"))
  
  var_sites <- map(new_samples, ~seg.sites(.x))
  
  samples_varsites <- map2(new_samples, var_sites, ~.x[,.y])
  
  var_sites_num <- map(var_sites, length) %>% 
    unlist()
  

  #cluster dictionary
  dictionary <- data.frame(unlist(file_list), var_sites_num) %>% 
    rename("file_name" = "unlist.file_list.")
  
  cluster_counts <- map2(snp_frame_level, id_times, ~.x %>%
                           dplyr::select(id, cluster, count) %>%
                           left_join(.y, by = c("id" = "SampleID")) %>%
                           rename("seq.name" = "date_id") %>%
                           mutate(cluster_name = str_c(snp_level,
                                                       "lineage", 
                                                       lineage_rough, 
                                                       "cluster", 
                                                       cluster, 
                                                       sep = "_")) %>%
                           ungroup() %>%
                           dplyr::select(cluster_name, count) %>%
                           unique()) %>% 
    bind_rows() %>% 
    mutate(file_name = str_c(cluster_name, ".fasta", sep=""))
  
  
  dictionary <- dictionary %>% 
    left_join(cluster_counts, by = "file_name") %>%
    mutate(snp_level = level)
  
  write_csv(dictionary, here::here("data", str_c("snp", level, "kopanyo_cluster_dictionary.csv", sep = "_")))
}

# splitting by SNP level five ---------------------------------------------

snp_frame_five <- map(snp_frame, ~.x %>%
                                 filter(snp_level == 5) %>%
                                 group_by(cluster) %>% 
                                 mutate(count = n()) %>% 
                                 filter(count >= 4))

samples_five <- map2(samples, snp_frame_five,  ~.x[rownames(.x) %in% .y$id,])




#let's do the whole whammy here, let's rename the sequences with their sample dates and then also grab 

samples1 <- read.fasta(file=here("data", "Botswana_lineage_1_86_samples.fasta"))

samples2 <- read.fasta(file=here("data", "Botswana_lineage_2_76_samples.fasta"))


samples3 <- read.fasta(file=here("data", "Botswana_lineage_3_12_samples.fasta"))

samples4 <- read.fasta(file=here("data", "Botswana_lineage_4_1252_samples.fasta"))


samples_data_frame <- list(samples1, samples2, samples3, samples4)


id_times <- map(meta_data, ~.x %>%
                            dplyr::select(SampleID, collectdt, lineage_rough) %>%
                            mutate(date_id = str_c(SampleID, as.character(collectdt), sep = "_")))



samples_data_frame <- pmap(list(samples_data_frame, id_times, snp_frame_five), ~..1 %>% 
                             left_join(..2, by = c("seq.name" = "SampleID")) %>%
                             filter(seq.name %in% ..3$id) %>%
                             dplyr::select(date_id, seq.text) %>%
                             rename("seq.name" = "date_id"))
                           
                           
                           
cluster_names <- map2(snp_frame_five, id_times, ~.x %>%
                        dplyr::select(id, cluster) %>%
                        left_join(.y, by = c("id" = "SampleID")) %>%
                        rename("seq.name" = "date_id") %>%
                        mutate(cluster_name = str_c("snp5", "lineage", lineage_rough, "cluster", cluster, sep = "_")) %>%
                        ungroup() %>%
                        dplyr::select(seq.name, cluster_name))

map2(samples_data_frame, cluster_names, split_dat)


# reducing down to varying sites ------------------------------------------
file_list <- list.files("C:/Users/fiddl/Documents/kopanyo-archived-phylo/data/", pattern = "^snp5")
new_samples <- map(file_list, ~read.dna(here::here("data", .x), format = "fasta"))

var_sites <- map(new_samples, ~seg.sites(.x))


test_samples <- new_samples[[1]]

test_var <- var_sites[[1]]

testing <- test_samples[,test_var]
samples_varsites <- map2(new_samples, var_sites, ~.x[,.y])

var_sites_num <- map(var_sites, length) %>% 
                 unlist()

quantile(var_sites_num)


cutoff <- var_sites_num[var_sites_num >8]

test <- read.dna(here::here("data", file_list[[39]]), format = "fasta")


final_samples <- samples_varsites[-18]
final_file_list <- file_list[-18]

write.dna()
map2(final_file_list, final_samples,
     ~write.dna(.y, here::here("data", str_c("varsites", .x, sep = "_")), format = "fasta"))

#cluster dictionary
dictionary <- data.frame(unlist(file_list), var_sites_num) %>% 
              rename("file_name" = "unlist.file_list.")

cluster_counts <- map2(snp_frame_five, id_times, ~.x %>%
                        dplyr::select(id, cluster, count) %>%
                        left_join(.y, by = c("id" = "SampleID")) %>%
                        rename("seq.name" = "date_id") %>%
                        mutate(cluster_name = str_c("snp5", "lineage", lineage_rough, "cluster", cluster, sep = "_")) %>%
                        ungroup() %>%
                        dplyr::select(cluster_name, count) %>%
                        unique()) %>% 
                   bind_rows() %>% 
                   mutate(file_name = str_c(cluster_name, ".fasta", sep=""))

dictionary <- dictionary %>% 
              left_join(cluster_counts, by = "file_name")

write_csv(dictionary, here::here("data", "kopanyo_cluster_dictionary.csv"))
