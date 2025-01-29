
# Simulate Original Data

# load libraries
library(tidyverse)



get_simdata <- function(n, beta_treat, error){
  beta0 <- 1
  x <- rbinom(n, 1, prob = 0.5)
  if (error == "normal") {
    epsilon <- rnorm(n, mean = 0, sd = 2)
  }
  else {
    epsilon <-rlnorm(n, meanlog = 0, sdlog = log(2))
  }
  
  y <- beta0 + beta_treat * x + epsilon
  
  tibble(
    x = x,
    y = y
  )
  
}
