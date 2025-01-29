
#Function to extract estimates

library(broom)
library(tictoc)

source(here::here("functions", "02_fit_model.R"))

get_beta_estimates <- function(model) {
    #if we have small enough sample size and all x values are the same, the model won't fit a treatment effect --> automatically 0
  ifelse(is.na(model$coefficients[[2]]), 0, model$coefficients[[2]])
}

get_waldCI_estimates <- function(model, beta_true, errors_true, alpha){
  tic() #start time
    
    #Wald Interval
  
  estimates <- tidy(model, conf.int = TRUE, conf.level = 1 - alpha) %>%
      filter(term == "x") %>%
      mutate(coverage = ifelse(beta_true >= conf.low & beta_true <= conf.high, 1, 0)) %>%
      rename(beta_hat = estimate, lower_CI = conf.low, upper_CI = conf.high) %>%
      select(beta_hat, lower_CI, upper_CI, coverage)
  
  time_end <- toc(quiet=TRUE) #end time
  
  estimates$time <- time_end$toc - time_end$tic #time elapsed
  
  return(estimates)
  
}

get_bootstrap_PI_estimates <- function(model, data, n, beta_true, errors_true, alpha, n_boot = 100){
  tic() #start time
  
    #Extract Beta Estimate from original sample
  beta_hat <- get_beta_estimates(model)

  
  #vector to save bootstrap samples in
  beta_hat_b <- rep(NA,n_boot) #beta
  
  
  #loop for outer bootstrap
  for (b in 1:n_boot) {
    
      #sample with replacement
    data_b <- data[sample(1:nrow(data), size = n, replace = TRUE),]
      #fit model to bootstrap sample
    model_b <- get_model_fit(data = data_b)
      #extract beta estimate
    beta_hat_b[[b]] <- get_beta_estimates(model_b)
  
  }
  
  #percentile interval
  per_int <- quantile(beta_hat_b, probs= c(alpha/2, 1-alpha/2))
  PI_coverage <- ifelse(per_int[[1]] < beta_true & per_int[[2]] > beta_true, 1, 0)
  
  estimates <- data.frame(beta_hat, per_int[[1]], per_int[[2]], PI_coverage)
  colnames(estimates) <- c("beta_hat", "lower_CI", "upper_CI", "coverage")
  estimates <- tibble(estimates)
  
  time_end <- toc(quiet=TRUE) #end time
  
  estimates$time <- time_end$toc - time_end$tic #time elapsed
  
  return(estimates)
  
}

get_bootstrap_t_estimates <- function(model, data, n, beta_true, errors_true, alpha, n_boot = 100, n_inner_boot = 50){
  tic() #start time
  
  #Extract Beta Estimate from original sample
  beta_hat <- get_beta_estimates(model)
  
  #vectors to save bootstrap samples in
  tstar_b <- rep(NA, n_boot) #t stat
  beta_hat_b <- rep(NA, n_boot) #beta
  
  
  #loop for outer bootstrap
  for (b in 1:n_boot) {
    
    #sample with replacement
    data_b <- data[sample(1:nrow(data), size = n, replace = TRUE),]
    #fit model to bootstrap sample
    model_b <- get_model_fit(data = data_b)
    #extract beta estimate
    beta_hat_b[[b]] <- get_beta_estimates(model_b)
    
    #create a vector to save inner bootstrap
    beta_hat_b_k <- rep(NA, n_inner_boot)
    #inner bootstrap
    for (k in 1:n_inner_boot) {
      
      #sample with replacement from outer bootstrap sample
      data_b_k <- data_b[sample(1:nrow(data_b), size = n, replace = TRUE),]
      #fit model to bootstrap sample
      model_b_k <- get_model_fit(data = data_b_k)
      #extract beta estimate
      beta_hat_b_k[[k]] <- get_beta_estimates(model_b_k)
    }
    
    #compute se as the sd of the inner bootstrap beta estimates
    se_beta_hat_b <- sd(beta_hat_b_k)
    #compute t for the outer bootstrap
    tstar_b[[b]] <- (beta_hat_b[[b]] - beta_true)/se_beta_hat_b
    
  }
  
  tstar_lower <- quantile(tstar_b, probs = 1 - alpha/2, na.rm = TRUE)[[1]]
  tstar_upper <- quantile(tstar_b, probs = alpha/2, na.rm = TRUE)[[1]]
  
  se_beta_hat <- sd(beta_hat_b)
  
  tI_lower <- beta_hat - tstar_lower*se_beta_hat
  tI_upper <- beta_hat - tstar_upper*se_beta_hat
  
  tI_coverage <- ifelse(tI_lower < beta_true & tI_upper > beta_true, 1, 0)
  
  estimates <- data.frame(beta_hat, tI_lower, tI_upper, tI_coverage)
  colnames(estimates) <- c("beta_hat", "lower_CI", "upper_CI", "coverage")
  estimates <- tibble(estimates)
  
  time_end <- toc(quiet=TRUE) #end time
  
  estimates$time <- time_end$toc - time_end$tic #time elapsed
  
  return(estimates)
  
}

