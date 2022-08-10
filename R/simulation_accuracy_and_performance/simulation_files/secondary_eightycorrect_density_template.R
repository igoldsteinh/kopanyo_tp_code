# Simulate TB-like outbreaks and generate timed phylogenies and TP posteriors
# Secondary simulation setting "increase sampling density" as described in manuscript
# Run for seeds 1 to 100 to replicate manuscript analysis
library(nosoi)  
library(tidyverse)
library(TransPhylo)
library(ape)
library(tidytree)
library(treeio)
library(rstan)
library(tidybayes)
library(coda)

args <- commandArgs(trailingOnly=TRUE)
print(args)
seed <- as.integer(args[1])
print(seed)

set.seed(as.numeric(seed))


# simulating epidemics ----------------------------------------------------

# pexit function
p_exit_fct  <- function(t, t_incub, p_ind_exit) {
  if(t < t_incub) {
    return(0)}
  if(t >= t_incub) {
    return(p_ind_exit)}
}

p_ind_exit_fct <- function(x) {
  rbeta(x, shape1 = 1/.11, shape2 = 3/.11)
}

# structure matrix moving between infector type
structure_matrix <- cbind(c(0,1), c(1,0))

colnames(structure_matrix) <- c("A", "B")

rownames(structure_matrix) <- c("A", "B")

# prob of being in each infector class
pmove_fct <- function(t, current.in) {
  if(t == 1) {
    if (current.in == "A") {
      return(.53)}
    if (current.in == "B") {
      return(.47)
    }
  }
  else {
    return(0)
  }
}

# probability of transmission given contact
p_trans_fct <- function(t, current.in, t_incub) {
  if (current.in == "A") {
    if (t < t_incub) {
      p = 0
    }
    if (t >= t_incub) {
      p = 0.0254467
    }
  }
  if (current.in == "B") {
    if (t<t_incub) {
      p = 0
    }
    if (t >= t_incub) {
      p = 0.01454097
    }
  }
  return(p)
}

# incubation time
t_incub_fct <- function(x) {
  
  t_incub <- rgamma(x, shape = 9 * .75, rate = .75)
  return(t_incub)
}

p_max_a_fct <- function(x) {
  rbeta(x, shape1 = 177.325, shape2 = 5611.62)
}

p_max_b_fct <- function(x) {
  rbeta(x, shape1 = 1.02/.05, shape2 = 99.02/.05) # 
}



# num contacts
n_contact_fct <- function(t) {
 
  contact <- abs(round(rnorm(1, 20, 5), 0))
  return(contact)
}

num_sims <- 10
sim_list <- vector(mode = "list", length = num_sims)
c <- 0
for (i in 1:num_sims) {
condition <- FALSE
while(condition == FALSE) {
test_sim <- nosoiSim(type = "single", 
                     popStructure = "discrete",
                     length.sim = 96, 
                     max.infected = 100000,
                     init.individuals = 1,
                     init.structure = "A",
                     structure.matrix = structure_matrix,
                     
                     pExit = p_exit_fct,
                     param.pExit = list(t_incub = t_incub_fct, 
                                        p_ind_exit = p_ind_exit_fct),
                     
                     pMove = pmove_fct,
                     param.pMove = NA,
                     timeDep.pMove = FALSE,
                     diff.pMove = TRUE,
                     
                     
                     nContact = n_contact_fct,
                     param.nContact = NA,
                     
                     pTrans = p_trans_fct,
                     param.pTrans = list(t_incub = t_incub_fct),
                     timeDep.pTrans = FALSE,
                     diff.pTrans = TRUE,
                     prefix.host = "H",
                     print.progress = FALSE)
if (test_sim[["host.info.A"]]$N.infected > 50 & test_sim[["host.info.A"]]$N.infected <= 2000) {
  condition <- TRUE
}
 c <- c + 1
}
sim_list[[i]] <- test_sim

}


# creating sample trees ---------------------------------------------------
tree_list <- map(sim_list, getTransmissionTree)

outputs <- map(sim_list, getTableHosts) %>%
           map(~.x %>% mutate(max_time = max(out.time[!is.na(out.time)]),
                              new_out = ifelse(is.na(out.time), max_time, out.time),
                              difftime = new_out - inf.time, 
                              sans_incub_time = difftime - t_incub)) %>%
            map(~.x %>% mutate(active_start = ceiling(inf.time + t_incub),
                                             active_time = new_out - active_start,
                                             last_one = new_out >= max_time - 11,
                                             last_two = new_out >= max_time - 23 & active_start < max_time - 11,
                                             last_three = new_out >= max_time - 35 & active_start < max_time - 23))

filtered_samp_one <- map(outputs, ~.x %>% filter(last_one == TRUE & active_time > 0) %>%
                           mutate(year_start = max_time - 11,
                                  year_end = max_time))

filtered_samp_two <- map(outputs, ~.x %>% filter(last_two == TRUE & active_time > 0) %>%
                           mutate(year_start = max_time - 23,
                                  year_end = max_time - 12 ))

filtered_samp_three <- map(outputs, ~.x %>% filter(last_three == TRUE & active_time > 0) %>%
                             mutate(year_start = max_time - 35,
                                    year_end = max_time - 24))

num_infectors_one <- map(filtered_samp_one, dim) %>%
  map(1)

num_infectors_two <- map(filtered_samp_two, dim) %>%
  map(1)


num_infectors_three <- map(filtered_samp_three, dim) %>%
  map(1)


num_samples_one <- map(num_infectors_one, ~round(.x*.8))
num_samples_two <- map(num_infectors_two, ~round(.x*.8))
num_samples_three <- map(num_infectors_three, ~round(.x*.8))

random_samp_one <- map2(filtered_samp_one, num_samples_one, sample_n)

filtered_samp_two <- map2(filtered_samp_two, random_samp_one, ~.x %>%
                            filter(!(hosts.ID %in% .y$hosts.ID)))


random_samp_two <- map2(filtered_samp_two, num_samples_two, sample_n)

filtered_samp_three <- map2(filtered_samp_three, random_samp_one, ~.x %>%
                              filter(!(hosts.ID %in% .y$hosts.ID))) %>%
  map2(random_samp_two, ~.x %>%
         filter(!(hosts.ID %in% .y$hosts.ID)))

random_samp_three <- map2(filtered_samp_three, num_samples_three, sample_n)


random_samp <- pmap(list(random_samp_one, random_samp_two, random_samp_three), rbind)


total_infectors <- map(outputs, dim)%>%
  map(1)

total_sampled <-  map(random_samp, dim)%>%
  map(1)
sampling_proportion <- map2(total_sampled,total_infectors,  ~.x/.y)
# write_rds(random_samp, here::here("code", "R Code", "hiv_cluster_random_samp_mean3_newsamp.rds"))

final_sample_table <- map(random_samp, ~.x %>% 
                            rowwise %>% 
                            mutate(sample_time = round(runif(1, active_start, new_out)),
                                   problems = sample_time > new_out | sample_time > out.time)%>%
                            dplyr::select(hosts.ID, sample_time) %>%
                            mutate(labels = paste("Sample", hosts.ID, sep = "_")) %>%
                            rename("hosts"= "hosts.ID", "times" = "sample_time"))

sample_inputs <- list(sim_list, tree_list, final_sample_table)




sample_tree <- pmap(sample_inputs, sampleTransmissionTree)

sample_tree <- map(sample_tree, as.phylo)

for (i in 1:length(sample_tree)) {
  sample_tree[[i]] <- multi2di(sample_tree[[i]])
  
  sample_tree[[i]]$edge.length <- pmax(sample_tree[[i]]$edge.length, 1/30)
  
  sample_tree[[i]]$edge.length <- sample_tree[[i]]$edge.length/12
}

date_last_sample <- map(outputs, ~.x %>% 
                          dplyr::select(max_time) %>% 
                          mutate(new_max_time = max_time/12) %>% 
                          unique %>% pull)
p_tree<-map2(sample_tree, date_last_sample, ptreeFromPhylo)

rm(tree_list)
rm(filtered_samp_one)
rm(filtered_samp_two)
rm(filtered_samp_three)
rm(random_samp_one)
rm(random_samp_two)
rm(random_samp_three)
rm(sample_inputs)
# run transphylo ----------------------------------------------------------

# generation time gamma distribution parameters
w.shape = 10
w.scale = 1/10

# sampling time gamma distribution parameters
ws.shape = 10
ws.scale = 1/10

# infer the transmission tree
res <- infer_multittree_share_param(p_tree,
                                    mcmcIterations = 100000,
                                    thinning = 10,
                                    w.shape = w.shape,
                                    w.scale = w.scale,
                                    ws.shape = ws.shape,
                                    ws.scale = ws.scale,
                                    updateOff.p = FALSE,
                                    startOff.p = .5,
                                    prior_pi_a = 1,
                                    prior_pi_b = 19,
                                    share = c("off.r", "neg", "off.p"),
                                    startNeg = 1.48,
                                    dateT = max(unlist(date_last_sample)) + .01)

rm(p_tree)
get_param_estimates <- function(record, p){
  sapply(record, function(x) x[[p]])
}


pi_param_results <- map(res, ~data.frame(pi = get_param_estimates(., "pi"))) %>%
  map(~.x %>% slice(5001:10000))



other_param_results <- map(res, ~data.frame(off.r = get_param_estimates(.x, "off.r"),
                                            off.p = get_param_estimates(.x, "off.p"),
                                            neg = get_param_estimates(.x, "neg"),
                                            pTTree = get_param_estimates(.x, "pTTree"),
                                            pPTree = get_param_estimates(.x, "pPTree")) %>%
                             mutate(log_posterior = pTTree + pPTree,
                                    R = off.r*off.p/(1-off.p)) %>%
                             slice(5001:10000))




coda_res <- map(res, convertToCoda)

coda_eff <- map(coda_res, effectiveSize) %>%
  bind_rows(.id = "column_label")


# models and prob inf distribution ----------------------------------------

names <- map(res, pluck, 1, "ctree", "nam")
offspring <- map2(res, names, getOffspringDist, burnin = .5) 
# times <-map2(res, names, getInfectionTimeDist, burnin = .5) %>%
  #      map(data.frame) %>%
   #     bind_rows(.id = "sim")

prob_inf <- map(offspring, ~1-rowSums(. == 0)/dim(.)[2])

prob_inf_frame <- map2(names, prob_inf, data.frame)%>%
  map(~.x %>% rename("names" =".x..i..", "prob_inf" = ".y..i.."))%>%
  map(~.x %>% mutate(new_names = str_replace(.$names, ".*_", "_"))) %>%
  map(~.x %>% mutate(final_names = str_replace_all(.$new_names, "_", "")))%>%
  map2(outputs, ~.x %>% mutate(true_infector = final_names %in% .y$inf.by,
                               tp_infector = prob_inf > .6)) 

hiv_status <- map(outputs, ~.x %>% dplyr::select(hosts.ID, current.in, inf.time, inf.by))

prob_inf_frame <- map2(prob_inf_frame, hiv_status, ~.x %>% 
                         left_join(.y, 
                                   by = c("final_names" = "hosts.ID"))) %>%
  map(~.x %>% dplyr::select(final_names, prob_inf, tp_infector, true_infector, current.in))


final_inf_frame <- bind_rows(prob_inf_frame, .id = "sim") %>%
  mutate(correct = tp_infector == true_infector)



model_tp <- glm(tp_infector ~ current.in, family = "binomial", data = final_inf_frame)

model_true <- glm(true_infector ~ current.in, family = binomial(link="logit"), data = final_inf_frame)




results <- list(model_tp,
                model_true,
                final_inf_frame,
                pi_param_results,
                other_param_results,
                coda_eff,
                sampling_proportion,
                sample_tree)

names(results) <- c("model_tp", 
                    "model_true", 
                    "final_inf_frame", 
                    "pi_param_results",
                    "other_param_results",
                    "coda_eff",
                    "sampling_proportion",
                    "sample_tree")

write_rds(results, str_c("secondary_eightycorrectdens_density_res_seed_", seed, ".rds", ""))

