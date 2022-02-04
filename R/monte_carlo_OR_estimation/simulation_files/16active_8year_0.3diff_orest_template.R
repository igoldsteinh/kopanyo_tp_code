# template for simming 1000 sims with pr transmission at 1.75 rather than 3
# Run for seeds 1-100
library(nosoi)  
library(tidyverse)
library(TransPhylo)
library(ape)
library(tidytree)
library(treeio)
library(rstan)
library(tidybayes)

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
      p = 0.008792846
      
    }
  }
  if (current.in == "B") {
    if (t<t_incub) {
      p = 0
    }
    if (t >= t_incub) {
      p = 0.02930949
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

num_sims <- 100
sim_list <- vector(mode = "list", length = num_sims)
c <- 0
for (i in 1:num_sims) {
  condition <- FALSE
  while(condition == FALSE) {
    test_sim <- nosoiSim(type = "single", 
                         popStructure = "discrete",
                         length.sim = 192, 
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
    if (test_sim[["host.info.A"]]$N.infected > 50) {
      condition <- TRUE
    }
    c <- c + 1
  }
  sim_list[[i]] <- test_sim
  
}

# tree_list <- map(sim_list, getTransmissionTree)

outputs <- map(sim_list, getTableHosts) %>%
  map(~.x %>% mutate(max_time = max(out.time[!is.na(out.time)]),
                     new_out = ifelse(is.na(out.time), max_time, out.time),
                     difftime = new_out - inf.time, 
                     sans_incub_time = difftime - t_incub))


outputs_table <- map(outputs, ~.x %>% mutate(active_start = ceiling(inf.time + t_incub),
                                             active_host = (!is.na(out.time) & active_start < out.time) | 
                                               (is.na(out.time) & active_start < max_time)))


odds <- map(outputs_table, ~.x %>%
              mutate(infector = hosts.ID %in% inf.by) %>%
              filter(active_host == TRUE) %>%
              group_by(current.in) %>%
              summarise(prob_inf = mean(infector),
                        prob_no_inf = 1-mean(infector),
                        num_total = n()) %>%
              ungroup()%>%
              mutate(odds_inf = prob_inf/prob_no_inf) %>%
              mutate(a_odds = lag(odds_inf),
                     odds_ratio = a_odds/odds_inf)
) %>%
  bind_rows(.id = "sim")

write_rds(odds, str_c("16active_8yearsim_0.3diff_odds_seed_", seed, ".rds", ""))
