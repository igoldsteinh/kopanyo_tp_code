# this is a file for functions which help process TP sim results
# and perform inference on TP labels
library(rstan)
library(tidyverse)
library(tidybayes)
library(SAMBA)


set.seed(1234)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# TP accuracy table function ----------------------------------------------

tp_acc <- function(res, trials, cutoff  = 0.6) {
  prob_inf <- map(res, pluck, "final_inf_frame") 

  prob_inf_frame <- prob_inf %>% 
    bind_rows(.id = "simsim") %>%
    mutate(tp_infector = prob_inf > cutoff,
           correct = tp_infector == true_infector)
  
  correctness <- prob_inf_frame %>%
    group_by(simsim, true_infector) %>%
    summarise(count_correct = sum(correct),
              total = n(), 
              percent = count_correct/total) %>%
    group_by(simsim) %>%
    mutate(num_cases = sum(total))
  
  mc_correctness <- correctness %>%
    group_by(true_infector) %>%
    summarise(
      mean_percent = mean(percent),
      se_percent = sd(percent)/sqrt(trials),
      mean_cases = mean(num_cases),
      sd_cases = sd(num_cases)
    ) 
  
  mc_correctness <- mc_correctness %>%
    mutate(
           num_trials = trials) %>% 
    dplyr::select(true_infector, mean_percent, se_percent, mean_cases, sd_cases, num_trials)
  
  mc_correctness
  
}


# final summary function --------------------------------------------------

final_summary <- function(res, 
                          spec, 
                          sens,  
                          trials, 
                          cluster = TRUE, 
                          true_val, 
                          shift = 0, 
                          add = TRUE,
                          stan = TRUE, 
                          cutoff = 0.6) {

  if (stan == TRUE){
    final_summary <- rbind(glm_summaries(res, true_val, cutoff), 
                           SAMBA_summary(res, true_val, shift, add, cutoff), 
                           stan_summary(res, spec, sens, cluster, true_val, shift, add, cutoff)) %>%
      mutate(num_trials = trials)
    
  }
  
  if (stan == FALSE) {
    final_summary <- rbind(glm_summaries(res, true_val, cutoff), 
                           SAMBA_summary(res, true_val, shift, add, cutoff)) %>%
      mutate(num_trials = trials)
    
  }
  
  return(final_summary)
}
# glm summaries -----------------------------------------------------------

glm_summaries <- function(res, true_val, cutoff) {
  if (cutoff == 0.6) {
    tp_glm <- map(res, pluck, "model_tp") %>%
      map(summary.glm)%>%
      map(~as.data.frame(.x[["coefficients"]]))%>%
      map(~.x %>% 
            tibble::rownames_to_column()%>%
            filter(rowname == "current.inB")) %>%
      bind_rows(.id = "sim") %>%
      mutate(odds_ratio = 1/exp(Estimate),
             below_one = odds_ratio < 1,
             reject = `Pr(>|z|)` < .05,
             est_low = Estimate - 1.96*`Std. Error`,
             est_high = Estimate + 1.96*`Std. Error`,
             or_high = 1/exp(est_low),
             or_low = 1/exp(est_high),
             reject_low = or_low > 1,
             reject_high = or_high < 1,
             coverage = true_val >=or_low & true_val <= or_high)
    
    
    tp_summary <- tp_glm %>%
      mutate(diff = abs(odds_ratio - true_val),
             bias = odds_ratio - true_val,
             sq_error = (odds_ratio - true_val)^2,
             ci_width = or_high - or_low) %>%
      summarise(
        mean_odds_ratio = round(mean(odds_ratio),2),
        mean_lb = round(mean(or_low), 2),
        mean_ub = round(mean(or_high), 2),
        mean_diff = mean(diff),
        mean_bias = mean(bias),
        sq_MSE = sqrt(mean(sq_error)),
        mean_ci_width = mean(ci_width),
        percent_coverage = mean(coverage),
        percent_reject = mean(reject),
        percent_reject_low = mean(reject_low),
        percent_reject_high = mean(reject_high)
      )%>%
      mutate(response = "transphylo")
    
    
    true_glm <- map(res, pluck, "model_true") %>%
      map(summary.glm)%>%
      map(~as.data.frame(.x[["coefficients"]]))%>%
      map(~.x %>% 
            tibble::rownames_to_column()%>%
            filter(rowname != "(Intercept)")) %>%
      bind_rows(.id = "sim") %>%
      mutate(odds_ratio = 1/exp(Estimate),
             below_one = odds_ratio < 1,
             reject = `Pr(>|z|)` < .05,
             est_low = Estimate - 1.96*`Std. Error`,
             est_high = Estimate + 1.96*`Std. Error`,
             or_high = 1/exp(est_low),
             or_low = 1/exp(est_high),           
             reject_low = or_low > 1,
             reject_high = or_high < 1,
             coverage = true_val>=or_low & true_val <= or_high)
    
    
    true_summary <- true_glm %>%
      mutate(diff = abs(odds_ratio - true_val),           
             bias = odds_ratio - true_val,
             sq_error = (odds_ratio - true_val)^2,
             ci_width = or_high - or_low) %>%
      summarise(
        mean_odds_ratio = round(mean(odds_ratio),2),
        mean_lb = round(mean(or_low), 2),
        mean_ub = round(mean(or_high), 2),
        mean_diff = mean(diff),      
        mean_bias = mean(bias),
        sq_MSE = sqrt(mean(sq_error)),
        mean_ci_width = mean(ci_width),
        percent_coverage = mean(coverage),
        percent_reject = mean(reject),
        percent_reject_low =mean(reject_low),
        percent_reject_high = mean(reject_high)
      )%>%
      mutate(response = "truth")
    
    
    num_samples <-  map(res, pluck, "final_inf_frame") %>%
      map(dim) %>%
      map(1) %>%
      unlist()
    
    full_summary2 <- rbind(tp_summary, true_summary) %>%
      mutate(mean_num_samples = mean(num_samples)) %>%
      dplyr::select(mean_OR = mean_odds_ratio, 
                    mean_lb,
                    mean_ub,
                    mean_diff, 
                    mean_bias,
                    sq_MSE,
                    percent_reject,
                    mean_ci_width,
                    percent_coverage,
                    percent_reject_low,
                    percent_reject_high,
                    type = response,
                    mean_num_samples)
    
    full_summary2
    
  }
  
  if (cutoff != 0.6) {
    prob_inf <- map(res, pluck, "final_inf_frame") %>%
      map(~.x %>%
            mutate(numeric_hiv = as.numeric(current.in == "A"),
                   tp_infector = prob_inf > cutoff))
    
    
    model_tp <- map(prob_inf, ~glm(tp_infector ~ current.in, family = "binomial", data = .x))
    
    model_true <- map(prob_inf, ~glm(true_infector ~ current.in, family = binomial(link="logit"), data = .x))
    
    
    tp_glm <- map(model_tp, summary.glm)%>%
      map(~as.data.frame(.x[["coefficients"]]))%>%
      map(~.x %>% 
            tibble::rownames_to_column()%>%
            filter(rowname == "current.inB")) %>%
      bind_rows(.id = "sim") %>%
      mutate(odds_ratio = 1/exp(Estimate),
             below_one = odds_ratio < 1,
             reject = `Pr(>|z|)` < .05,
             est_low = Estimate - 1.96*`Std. Error`,
             est_high = Estimate + 1.96*`Std. Error`,
             or_high = 1/exp(est_low),
             or_low = 1/exp(est_high),
             reject_low = or_low > 1,
             reject_high = or_high < 1,
             coverage = true_val >=or_low & true_val <= or_high)
    
    
    tp_summary <- tp_glm %>%
      mutate(diff = abs(odds_ratio - true_val),
             bias = odds_ratio - true_val,
             sq_error = (odds_ratio - true_val)^2,
             ci_width = or_high - or_low) %>%
      summarise(
        mean_odds_ratio = round(mean(odds_ratio),2),
        mean_lb = round(mean(or_low), 2),
        mean_ub = round(mean(or_high), 2),
        mean_diff = mean(diff),
        mean_bias = mean(bias),
        sq_MSE = sqrt(mean(sq_error)),
        mean_ci_width = mean(ci_width),
        percent_coverage = mean(coverage),
        percent_reject = mean(reject),
        percent_reject_low = mean(reject_low),
        percent_reject_high = mean(reject_high)
      )%>%
      mutate(response = "transphylo")
    
    
    true_glm <- map(model_true, summary.glm)%>%
      map(~as.data.frame(.x[["coefficients"]]))%>%
      map(~.x %>% 
            tibble::rownames_to_column()%>%
            filter(rowname != "(Intercept)")) %>%
      bind_rows(.id = "sim") %>%
      mutate(odds_ratio = 1/exp(Estimate),
             below_one = odds_ratio < 1,
             reject = `Pr(>|z|)` < .05,
             est_low = Estimate - 1.96*`Std. Error`,
             est_high = Estimate + 1.96*`Std. Error`,
             or_high = 1/exp(est_low),
             or_low = 1/exp(est_high),           
             reject_low = or_low > 1,
             reject_high = or_high < 1,
             coverage = true_val>=or_low & true_val <= or_high)
    
    
    true_summary <- true_glm %>%
      mutate(diff = abs(odds_ratio - true_val),           
             bias = odds_ratio - true_val,
             sq_error = (odds_ratio - true_val)^2,
             ci_width = or_high - or_low) %>%
      summarise(
        mean_odds_ratio = round(mean(odds_ratio),2),
        mean_lb = round(mean(or_low), 2),
        mean_ub = round(mean(or_high), 2),
        mean_diff = mean(diff),      
        mean_bias = mean(bias),
        sq_MSE = sqrt(mean(sq_error)),
        mean_ci_width = mean(ci_width),
        percent_coverage = mean(coverage),
        percent_reject = mean(reject),
        percent_reject_low =mean(reject_low),
        percent_reject_high = mean(reject_high)
      )%>%
      mutate(response = "truth")
    
    
    num_samples <-  map(res, pluck, "final_inf_frame") %>%
      map(dim) %>%
      map(1) %>%
      unlist()
    
    full_summary2 <- rbind(tp_summary, true_summary) %>%
      mutate(mean_num_samples = mean(num_samples)) %>%
      dplyr::select(mean_OR = mean_odds_ratio, 
                    mean_lb,
                    mean_ub,
                    mean_diff, 
                    mean_bias,
                    sq_MSE,
                    percent_reject,
                    mean_ci_width,
                    percent_coverage,
                    percent_reject_low,
                    percent_reject_high,
                    type = response,
                    mean_num_samples)
    
    full_summary2
    
    
    
  }
  
}


# SAMBA processing --------------------------------------------------------
# note this will use the SAMBA calculated sensitivity score
SAMBA_summary <- function(res, true_val, shift, add, cutoff) {
  prob_inf <- map(res, pluck, "final_inf_frame") %>%
    map(~.x %>%
          mutate(numeric_hiv = as.numeric(current.in == "A"),
                 tp_infector = prob_inf > cutoff))
  
  if (add == TRUE){
  SAMBA_sens <- map(prob_inf, ~sensitivity(as.numeric(.x$tp_infector), 
                                                     X = .x$numeric_hiv, 
                                                     mean(.x$true_infector) + shift, 
                                                     r = NULL,
                                                     weights = NULL)) %>%
    map(pluck(1))
  }
  
  if (add == FALSE) {
    SAMBA_sens <- map(prob_inf, ~sensitivity(as.numeric(.x$tp_infector), 
                                             X = .x$numeric_hiv, 
                                             mean(.x$true_infector) - shift, 
                                             r = NULL,
                                             weights = NULL)) %>%
      map(pluck(1))
    
  }

  SAMBA_res <- map2(prob_inf, SAMBA_sens, ~approxdist(as.numeric(.x$tp_infector),
                                                                .x$numeric_hiv, .y))
  
  SAMBA_res <- SAMBA_res %>%
    bind_rows(.id = "sim") %>%
    mutate(CI_Low = param - 1.96 * sqrt(variance),
           CI_High = param + 1.96 * sqrt(variance),
           OR = exp(param),
           OR_Low = exp(CI_Low),
           OR_High = exp(CI_High),
           reject_low = OR_Low > 1,
           reject_high = OR_High < 1,
           reject = !(OR_Low < 1 & OR_High > 1),
           Coverage = true_val >= OR_Low & true_val <= OR_High)
  
  SAMBA_res <- SAMBA_res %>%
    mutate(CI_width = OR_High - OR_Low,
           diff = abs(OR - true_val),
           bias = OR - true_val,
           sq_error = (OR - true_val)^2
    )
  
  
  num_samples <-  map(res, pluck, "final_inf_frame") %>%
    map(dim) %>%
    map(1) %>%
    unlist()
  
  
  SAMBA_res_summary <- SAMBA_res %>%
    summarise(
      mean_OR = mean(OR),
      mean_lb = mean(OR_Low),
      mean_ub = mean(OR_High),
      mean_diff = mean(diff),
      mean_bias = mean(bias),
      sq_MSE = sqrt(mean(sq_error)),
      percent_reject = mean(reject),
      mean_ci_width = mean(CI_width),
      percent_reject_low = mean(reject_low),
      percent_reject_high = mean(reject_high),
      percent_coverage = mean(Coverage),
      type = "SAMBA Calc",
      mean_num_samples = mean(num_samples)
    )
  
  SAMBA_res_summary

}

# STAN homebrew analysis --------------------------------------------------

stan_summary <- function(res, spec, sens, cluster = TRUE, true_val, shift, add, cutoff) {
  prob_inf <- map(res, pluck, "final_inf_frame") %>%
    map(~.x %>%
          mutate(numeric_hiv = as.numeric(current.in == "A"),
                 tp_infector = prob_inf > cutoff))
  
  if (add == TRUE){
    model_objects_test <-map(prob_inf, ~list(N = dim(.x)[1],
                                             z = .x$tp_infector,
                                             x = .x$numeric_hiv,
                                             specificity = spec + shift,
                                             sensitivity = sens + shift))
    
  }
  
  if (add == FALSE){
    model_objects_test <-map(prob_inf, ~list(N = dim(.x)[1],
                                             z = .x$tp_infector,
                                             x = .x$numeric_hiv,
                                             specificity = spec - shift,
                                             sensitivity = sens - shift))
    
    
  }
  
  
  if (cluster == FALSE) {
    test_fit <- map(model_objects_test, ~stan(file =here::here("R",
                                                               "misclass_logistic_regression.stan"),
                                              data = .x,
                                              seed = 45,
                                              iter = 2000,
                                              chains = 4))
    
  }
  
  else {
    test_fit <- map(model_objects_test, ~stan(file ="misclass_logistic_regression.stan",
                                              data = .x,
                                              seed = 45,
                                              iter = 2000,
                                              chains = 4))
  }
  
  draws <- map(test_fit, tidy_draws) %>%
    bind_rows(.id = "sim") %>%
    dplyr::select(sim, beta0, beta1) %>%
    group_by(sim) %>%
    median_qi()
  
  misclass_results <- draws %>%
    mutate(OR = exp(beta1),
           OR_Low = exp(beta1.lower),
           OR_High = exp(beta1.upper))
  
  
  
  num_samples <-  map(res, pluck, "final_inf_frame") %>%
    map(dim) %>%
    map(1) %>%
    unlist()
  
  stan_res_summary <- misclass_results %>%
    mutate(Coverage = true_val >= OR_Low & true_val <= OR_High,
           reject_low = OR_Low > 1,
           reject_high = OR_High < 1,
           reject = !(OR_Low < 1 & 1 < OR_High),
           bias = OR - true_val, 
           sq_error = (OR - true_val)^2,
           CI_width = OR_High - OR_Low,
           diff = abs(OR - true_val)) %>%
    summarise(
      mean_OR = mean(OR),
      mean_lb = mean(OR_Low),
      mean_ub = mean(OR_High),
      mean_diff = mean(diff),
      mean_bias = mean(bias),
      sq_MSE = sqrt(mean(sq_error)),
      percent_reject = mean(reject),
      mean_ci_width = mean(CI_width),
      percent_reject_low = mean(reject_low),
      percent_reject_high = mean(reject_high),
      percent_coverage = mean(Coverage),
      type = "Stan",
      mean_num_samples = mean(num_samples)
    )
  
  stan_res_summary
  
}


# delta calculations ------------------------------------------------------

delta_calc <- function(prob_a, prob_b, sd_a, sd_b, n) {
  cov_matrix <- matrix(1:4, nrow = 2, ncol = 2)
  cov_matrix[1,1] <- sd_a^2
  cov_matrix[1,2] <- 0
  cov_matrix[2,1] <- 0
  cov_matrix[2,2] <- sd_b^2
  
  grad <- c(0,0)
  x <- prob_a
  y <-prob_b
  grad[1] <- (1/(1-x)^2) * ((1-y)/y)
  grad[2] <- (x/(1-x)) * (-1/y^2)
  
  sd <- sqrt(t(grad)%*%cov_matrix%*%grad)
  se <- sd/sqrt(n)
  
  estimate <- (x/(1-x))/(y/(1-y))
  
  results <- c(estimate, se, n)
  
  return(results)
  
}


# glm tables -----------------------------------------------------------

glm_tables <- function(res, true_val) {
  tp_glm <- map(res, pluck, "model_tp") %>%
    map(summary.glm)%>%
    map(~as.data.frame(.x[["coefficients"]]))%>%
    map(~.x %>% 
          tibble::rownames_to_column()%>%
          filter(rowname == "current.inB")) %>%
    bind_rows(.id = "sim") %>%
    mutate(odds_ratio = 1/exp(Estimate),
           below_one = odds_ratio < 1,
           reject = `Pr(>|z|)` < .05,
           est_low = Estimate - 1.96*`Std. Error`,
           est_high = Estimate + 1.96*`Std. Error`,
           est_low80 = Estimate - 1.28*`Std. Error`,
           est_high80 = Estimate + 1.28*`Std. Error`,
           or_high = 1/exp(est_low),
           or_low = 1/exp(est_high),
           or_high80 = 1/exp(est_low80),
           or_low80 = 1/exp(est_high80),
           reject_low = or_low > 1,
           reject_high = or_high < 1,
           coverage = true_val >=or_low & true_val <= or_high,
           coverage80 = true_val >= or_low80 & true_val <= or_high80,
           model = "TP")
  
  
  true_glm <- map(res, pluck, "model_true") %>%
    map(summary.glm)%>%
    map(~as.data.frame(.x[["coefficients"]]))%>%
    map(~.x %>% 
          tibble::rownames_to_column()%>%
          filter(rowname != "(Intercept)")) %>%
    bind_rows(.id = "sim") %>%
    mutate(odds_ratio = 1/exp(Estimate),
           below_one = odds_ratio < 1,
           reject = `Pr(>|z|)` < .05,
           est_low = Estimate - 1.96*`Std. Error`,
           est_high = Estimate + 1.96*`Std. Error`,           
           est_low80 = Estimate - 1.28*`Std. Error`,
           est_high80 = Estimate + 1.28*`Std. Error`,
           or_high = 1/exp(est_low),
           or_low = 1/exp(est_high),
           or_high80 = 1/exp(est_low80),
           or_low80 = 1/exp(est_high80),
           reject_low = or_low > 1,
           reject_high = or_high < 1,
           coverage = true_val>=or_low & true_val <= or_high,
           coverage80 = true_val >= or_low80 & true_val <= or_high80,
           model = "True")
  
  

  num_samples <-  map(res, pluck, "final_inf_frame") %>%
    map(dim) %>%
    map(1) %>%
    unlist()
  
  full_summary2 <- rbind(tp_glm, true_glm) 
  
  full_summary2
  
}


# 2x2 counts for sanghyuk -------------------------------------------------

count_2x2 <- function(prob_inf_list) {
  count_labels <-map(prob_inf_list, ~.x %>%
    group_by(current.in) %>%
    summarise(num_total = n(),
              num_trueIS = sum(true_infector),
              num_trueNS = sum(true_infector== FALSE),
              num_tpIS = sum(tp_infector),
              num_tpNS = sum(tp_infector == FALSE),
              true_prob_IS = num_trueIS/num_total
    )) %>%
    bind_rows(.id = "sim")
  
  return(count_labels)
}

# true glm tables -----------------------------------------------------------

true_glm_tables <- function(res, true_val) {

  true_glm <- map(res, pluck, "model_true") %>%
    map(summary.glm)%>%
    map(~as.data.frame(.x[["coefficients"]]))%>%
    map(~.x %>% 
          tibble::rownames_to_column()%>%
          filter(rowname != "(Intercept)")) %>%
    bind_rows(.id = "sim") %>%
    mutate(odds_ratio = 1/exp(Estimate),
           below_one = odds_ratio < 1,
           reject = `Pr(>|z|)` < .05,
           est_low = Estimate - 1.96*`Std. Error`,
           est_high = Estimate + 1.96*`Std. Error`,           
           est_low80 = Estimate - 1.28*`Std. Error`,
           est_high80 = Estimate + 1.28*`Std. Error`,
           or_high = 1/exp(est_low),
           or_low = 1/exp(est_high),
           or_high80 = 1/exp(est_low80),
           or_low80 = 1/exp(est_high80),
           reject_low = or_low > 1,
           reject_high = or_high < 1,
           coverage = true_val>=or_low & true_val <= or_high,
           coverage80 = true_val >= or_low80 & true_val <= or_high80,
           model = "True")
  
  
  
  num_samples <-  map(res, pluck, "samples") %>%
    map(~.x %>% ungroup() %>% summarise(num_samples = n())) %>%
    bind_rows(.id = "sim")
  
  
  true_glm <- true_glm %>%
              left_join(num_samples, by = "sim")
  
  true_glm
}

# all probs from 2x2 table -------------------------------------------------

calc_all_probs <- function(prob_inf_list) {
  total_probs <-map(prob_inf_list, ~.x %>%
                       summarise(num_total = n(),
                                 num_trueIS = sum(true_infector),
                                 num_HIV = sum(current.in == "B"),
                                 true_prob_IS = num_trueIS/num_total,
                                 prob_HIV = num_HIV/ num_total
                       )) %>%
    bind_rows(.id = "sim") %>%
    ungroup() %>%
    dplyr::select(sim, true_prob_IS, prob_HIV)
  
  
  cond_IS <- map(prob_inf_list, ~.x %>%
                   group_by(current.in) %>%
                   summarise(num_total = n(),
                             num_trueIS = sum(true_infector),
                             cond_prob_IS = num_trueIS/num_total
                   )) %>%
    bind_rows(.id = "sim") %>%
    ungroup() %>%
    dplyr::select(sim, current.in, cond_prob_IS) %>%
    pivot_wider(names_from = current.in, values_from = cond_prob_IS) %>%
    rename("IS_HIVneg" = "A", "IS_HIVpos" = "B")
  
  cond_HIV <- map(prob_inf_list, ~.x %>%
                    group_by(true_infector) %>%
                    summarise(num_total = n(),
                              num_HIV = sum(current.in == "B"),
                              cond_prob_HIV = num_HIV/num_total
                    )) %>%
    bind_rows(.id = "sim") %>%
    ungroup() %>%
    dplyr::select(sim, true_infector, cond_prob_HIV) %>%
    pivot_wider(names_from = true_infector, values_from = cond_prob_HIV) %>%
    rename("HIV_ISneg" = "FALSE", "HIV_ISpos" = "TRUE")
  
  
  probs_frame <- total_probs %>%
                 left_join(cond_IS, by = "sim") %>%
                 left_join(cond_HIV, by = "sim")
  
  
  return(probs_frame)
}
