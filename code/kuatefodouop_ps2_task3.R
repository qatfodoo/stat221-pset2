source("poissonLogN_MCMC.R")
source("kuatefodouop_ps2_task2.R")

## Task 3: Evaluate coverage for a simple case

# Constants
J <- 1000
N <- 2
w <- rep(1, J) # equal weights

# Simulation parameters
B <-200 #1000 # Number of simulations
B.theta <- 20 #40 # theta drawss
B.y <- floor(B / B.theta) # y draws for each theta

mu.array <- c(1.6, 2.5, 5.2, 4.9)
sigsq.array <- c(0.7^2, 1.3^2, 1,3^2, 1.6^2)

# Run simulation

if (Sys.getenv("SLURM_JOB_ID") != "") { # Divide computation per tasks
  
  job.id <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
  print(paste("Job id", job.id, sep=": "))
  task.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  print(paste("Task id", task.id, sep=": "))
  param.id <- task.id %% 4 + 1 # Parameter handled by task
  
  gc() # garbage collection
  
  # log.theta and corresponding coverage
  log.theta_mat <- matrix(data=0, nrow=J, ncol= B.theta) # Store log.theta values
  log.theta_mean <- replicate(B.theta, matrix(data=0, nrow=J, ncol= B.y), simplify=F) # Store posterior means of draws
  log.theta_sd <- replicate(B.theta, matrix(data=0, nrow=J, ncol= B.y), simplify=F) # Store posterior sd of drawslog.theta_mean
  cov68_mat <- matrix(data=0, nrow=J, ncol= B.theta) # Store 68% frequency coverage for thetas
  cov95_mat <- matrix(data=0, nrow=J, ncol= B.theta) # Store 95% frequency coverage for theta
  
  # Set current model constant hyperparameters
  mu <- mu.array[param.id]
  sigsq <- sigsq.array[param.id]
  
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
    dt.mcmc <- (t2.mcmc - t1.mcmc) / 60 # dt in min
    print(paste("Current theta MCMCs elapsed time (min)", dt.mcmc, sep=": "))
    
    # Compute frequency coverage for log.theta
    cov68_mat[, b.theta] <- apply(logTheta.cov_68, 1, mean)
    cov95_mat[, b.theta] <- apply(logTheta.cov_68, 1, mean)
    
  }
  t2.sim <- Sys.time()
  dt.sim <- (t2.sim - t1.sim) / 60 # dt in min
  print(paste(paste("Simulation elapsed time (min), task", task.id, sep=" "), dt.sim, sep=": "))
  
  gc() # required memory
  
  # Store results in output folder
  save(list=c("log.theta_mat", "log.theta_mean", "log.theta_sd",
              "cov68_mat", "cov95_mat"), file=paste("./out/task3_out_jobid", job.id, "_taskid",
                task.id, "_param", param.id, ".Rdata", sep=""))
}


