source("poissonLogN_MCMC.R")
source("kuatefodouop_ps2_task2.R")

## Task 3: Evaluate coverage for a simple case

# Constants
J <- 1000
N <- 2
w <- rep(1, J) # equal weights

# Simulation parameters
B <-1 #1000 # Number of simulations
B.theta <- 1 #40 # theta drawss
B.y <- floor(B / B.theta) # y draws for each theta

mu.array <- c(1.6, 2.5, 5.2, 4.9)
sigsq.array <- c(0.7^2, 1.3^2, 1,3^2, 1.6^2)

# log.theta and corresponding coverage
log.theta_mat <- matrix(data=0, nrow=J, ncol= B.theta) # Store log.theta values
log.theta_mean <- replicate(B.theta, matrix(data=0, nrow=J, ncol= B.y), simplify=F) # Store posterior means of draws
log.theta_sd <- replicate(B.theta, matrix(data=0, nrow=J, ncol= B.y), simplify=F) # Store posterior sd of drawslog.theta_mean
cov68_mat <- matrix(data=0, nrow=J, ncol= B.theta) # Store 68% frequency coverage for thetas
cov95_mat <- matrix(data=0, nrow=J, ncol= B.theta) # Store 95% frequency coverage for theta

t1.ov <- Sys.time()
for (i in 1:1) { #length(mu)) {
  
  # Set current model constant hyperparameters
  mu <- mu.array[i]
  sigsq <- sigsq.array[i]
  
  # Draw set of thetas
  t1.sim <- Sys.time()
  for (b.theta in 1:B.theta) {
    
    log.theta <- rnorm(J, mu, sqrt(sigsq)) # Sample log.theta
    log.theta_mat[, b.theta] <- log.theta
    theta <- exp(log.theta)
    
    logTheta.cov_68 <- matrix(data=0, nrow=J, ncol= B.y) # Store 68% coverage for thetas
    logTheta.cov_95 <- matrix(data=0, nrow=J, ncol= B.y) # Store 95% coverage for theta
    
    
    # Draw set of Y for current theta
    t1.mcmc <- Sys.time()
    for (b.y in 1:B.y) {
      
      Y <- simYgivenTheta(theta, w, N) # Simulate observations
      mcmc <- poisson.logn.mcmc(Y, w)
      
      logTheta.draws <- mcmc$logTheta # Sample draws
      log.theta_mean[[b.theta]][, b.y] <- apply(logTheta.draws, 1, mean) # Store sample means
      log.theta_sd[[b.theta]][, b.y] <- apply(logTheta.draws, 1, sd) # Store sample sd
      
      # Derive 68% and 95% poterior intervals
      post.int_68 <- t(apply(logTheta.draws, 1, function (r) { quantile(r, probs=c(0.16, 0.84)) }))
      post.int_95 <- t(apply(logTheta.draws, 1, function (r) { quantile(r, probs=c(0.025, 0.975)) }))
      
      # Store current coverage
      logTheta.cov_68[, b.y] <- (log.theta >= post.int_68[, 1]) & (log.theta <= post.int_68[, 2])
      logTheta.cov_95[, b.y] <- (log.theta >= post.int_95[, 1]) & (log.theta <= post.int_95[, 2])
    }
    t2.mcmc <- Sys.time()
    dt.mcmc <- t2.mcmc - t1.mcmc
    print(paste("Current theta MCMCs elapsed time", dt.mcmc, sep=": "))
    
    # Compute frequency coverage for log.theta
    cov68_mat[, b.theta] <- apply(logTheta.cov_68, 1, mean)
    cov95_mat[, b.theta] <- apply(logTheta.cov_68, 1, mean)
    
  }
  t2.sim <- Sys.time()
  dt.sim <- t2.sim - t1.sim
  print(paste("Simulation elapsed time", dt.sim, sep=": "))
}
t2.ov <- Sys.time()
dt.ov <- t2.ov - t1.ov
print(paste("Total elapsed time", dt.ov, sep=": "))
