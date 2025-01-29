
  #Run Simulation from scenario params

get_simulated_scenario <- function(params, alpha) {
    #Set Seeds for looping
  seed <- floor(runif(params$n_sim, 1, 10000))
  
  # Set up Parallel back end with 6 cores
  num_cores <- detectCores() - 2
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  parallel_results <- foreach(i = 1:params$n_sim, .packages = c("tidyverse", "broom", "tictoc"), .combine = rbind) %dopar% {
    set.seed(seed[i])
    
    source(here::here("functions", "01_simulate_data.R"))
    source(here::here("functions", "02_fit_model.R"))
    source(here::here("functions", "03_extract_estimates.R"))
    
    ####################
    # simulate original data
    sim_data <- get_simdata(n = params$n,
                            beta_treat = params$beta_true,
                            error = params$errors_true)
    
    ####################
    # apply method(s)
    sim_model <- get_model_fit(data = sim_data)
    
    
    ####################
    # calculate estimates
    sim_waldCI <- get_waldCI_estimates(model = sim_model,
                                       beta_true = params$beta_true,
                                       errors_true = params$errors_true, 
                                       alpha = alpha)
    
    sim_bootstrap_PI <- get_bootstrap_PI_estimates(model = sim_model,
                                                   data = sim_data,
                                                   n = params$n,
                                                   beta_true = params$beta_true,
                                                   errors_true = params$errors_true,
                                                   alpha = alpha, 
                                                   n_boot = params$n_boot)
    
    sim_bootstrap_tI <- get_bootstrap_t_estimates(model = sim_model,
                                                  data = sim_data,
                                                  n = params$n,
                                                  beta_true = params$beta_true,
                                                  errors_true = params$errors_true,
                                                  alpha = alpha, 
                                                  n_boot = params$n_boot,
                                                  n_inner_boot = params$n_inner_boot)
    
    ####################
    # store results, including estimates, speed, parameter scenarios
    result <- data.frame(params,
                         seed = seed[i], 
                         beta_hat = sim_waldCI$beta_hat, 
                         Wald_lower = sim_waldCI$lower_CI,
                         Wald_upper = sim_waldCI$upper_CI,
                         Wald_coverage = sim_waldCI$coverage,
                         Wald_time = sim_waldCI$time,
                         BS_PI_lower = sim_bootstrap_PI$lower_CI,
                         BS_PI_upper = sim_bootstrap_PI$upper_CI,
                         BS_PI_coverage = sim_bootstrap_PI$coverage,
                         BS_PI_time = sim_bootstrap_PI$time,
                         BS_tI_lower = sim_bootstrap_tI$lower_CI,
                         BS_tI_upper = sim_bootstrap_tI$upper_CI,
                         BS_tI_coverage = sim_bootstrap_tI$coverage,
                         BS_tI_time = sim_bootstrap_tI$time)
    
    result
    
  }
  
  return(parallel_results)
  
}

