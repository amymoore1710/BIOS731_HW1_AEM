
#Function to Fit Model

get_model_fit <- function(data){
  lm(y ~ x, data = data)
}
