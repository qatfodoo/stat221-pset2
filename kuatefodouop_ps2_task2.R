## Task 2: Write functions to simulate data from the model

simYgivenTheta <- function(theta, w, N) {
  # Parameters of the j poisson dist
  lambda = w * theta
  
  t(sapply(lambda, function(x) {rpois(N, x)}))
}
