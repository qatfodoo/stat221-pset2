source("poissonLogN_MCMC.R")
source("kuatefodouop_ps2_task2.R")

## Task 3: Evaluate coverage for a simple case

# Constants
J <- 1000
N <- 2
weights <- rep(1, J) # equal weights

# Simulation parameters
B <- 1000 # Number of simulations
B.theta <- 40 # theta drawss
B.y <- 25 # y draws for each theta

mu <- c(1.6, 2.5, 5.2, 4.9)
sigsq <- c(0.7^2, 1.3^2, 1,3^2, 1.6^2)