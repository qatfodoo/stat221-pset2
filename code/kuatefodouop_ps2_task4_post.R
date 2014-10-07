library(ggplot2)
source("multiplot.R")

## Task 4 post-processing

# Simulation parameters for object initialization
J <- 1000
N <- 2
w <- as.numeric(unlist(read.table("weights.txt"))) # distinct weights

B <-40 # Number of simulations
B.theta <- 4 # theta drawss
B.y <- floor(B / B.theta) # y draws for each theta
n.param <- 4 # Number of parameters
n.t <- 3 # Number of tasks per parameter

# Initialize list of parameter's stacked outputs
ltheta.list <- replicate(n.param, matrix(data=0, nrow=J, ncol= B.theta * n.t), simplify=F) # Store log.theta values
ltheta_mean.list <- replicate(n.param, replicate(B.theta * n.t, matrix(data=0, nrow=J, ncol= B.y),
                                               simplify=F), simplify=F) # Store posterior means of draws
ltheta_sd.list <- replicate(n.param, replicate(B.theta * n.t, matrix(data=0, nrow=J, ncol= B.y),
                                               simplify=F), simplify=F) # Store posterior sd of drawslog.theta_mean
cov68.list <- replicate(n.param, matrix(data=0, nrow=J, ncol= B.theta * n.t), simplify=F) # Store 68% frequency coverage for thetas
cov95.list <- replicate(n.param, matrix(data=0, nrow=J, ncol= B.theta * n.t), simplify=F) # Store 95% frequency coverage for theta
logw.list <- replicate(n.param, matrix(data=log(w), nrow=J, ncol= B.theta * n.t), simplify=F)

# Get list of files corresponding to each parameter
files.param <- sapply(1:4, function(i) {
                            list.files(path = "./out/task4", pattern = paste("param" ,i, "+", sep=""))
                              }, simplify=F)

# Assemble output files from different tasks

for (p in 1:n.param) {
  
  for (t in 1:3) {
    
    # Load output from current task for parameter
    load(paste("./out/task4/", files.param[[p]][t], sep=""))
    
    # Update stacked output
    ltheta.list[[p]][, ((t - 1) * B.theta + 1):(t * B.theta)] <- log.theta_mat
    for (i.logt in 1:B.theta) {
      ltheta_mean.list[[p]][[(t - 1) * B.theta + i.logt]] <- log.theta_mean[[i.logt]]
      ltheta_sd.list[[p]][[(t - 1) * B.theta + i.logt]] <- log.theta_sd[[i.logt]]
    }
    cov68.list[[p]][, ((t - 1) * B.theta + 1):(t * B.theta)] <- cov68_mat
    cov95.list[[p]][, ((t - 1) * B.theta + 1):(t * B.theta)] <- cov95_mat
    
  }
}

# Plot coverage against parameters

p_logt68.list <- vector("list", n.param) # List of plots against log theta
p_logt95.list <- vector("list", n.param) # List of plots against log theta
p_logw68.list <- vector("list", n.param) # List of plots against log w
p_logw95.list <- vector("list", n.param) # List of plots against log w

for (p in 1:n.param) {
  
  ltheta <- ltheta.list[[p]]
  ltheta_mean <- ltheta_mean.list[[p]]
  ltheta_sd <- ltheta_sd.list[[p]]
  cov68 <- cov68.list[[p]] 
  cov95 <- cov95.list[[p]]
  logw <- logw.list[[p]]
  
  # Plot against log theta
  x11() # Create new plot window

  p_logt68.list[[p]] <- ggplot(, aes(x=as.vector(ltheta), as.vector(cov68))) + 
    geom_point(size=3, alpha=1, colour="darkblue") +
    geom_smooth(aes(group=1), method = 'loess', size=2, colour='red', se=FALSE) +
    ylim(0, 1) +
    ylab("cov_68") +
    xlab("log.theta") +
    labs(title="68% coverage against Log Theta, task 4")
  
  p_logt95.list[[p]] <- ggplot(, aes(x=as.vector(ltheta), as.vector(cov95))) + 
    geom_point(size=3, alpha=1, colour="darkblue") +
    geom_smooth(aes(group=1), method = 'loess', size=2, colour='red', se=FALSE) +
    ylim(0, 1) +
    ylab("cov_95") +
    xlab("log.theta") +
    labs(title="95% coverage against Log Theta, task 4")
  
  multiplot(p_logt68.list[[p]], p_logt95.list[[p]])
  
  # Plot againt log w
  x11() # Create new plot window
  
  p_logw68.list[[p]] <- ggplot(, aes(x=as.vector(logw), as.vector(cov68))) + 
    geom_point(size=3, alpha=1, colour="darkblue") +
    geom_smooth(aes(group=1), method = 'loess', size=2, colour='red', se=FALSE) +
    ylim(0, 1) +
  ylab("cov_68") +
    xlab("log.w") +
    labs(title="68% coverage against Log w, task 4")
  
  p_logw95.list[[p]] <- ggplot(, aes(x=as.vector(logw), as.vector(cov68))) + 
    geom_point(size=3, alpha=1, colour="darkblue") +
    geom_smooth(aes(group=1), method = 'loess', size=2, colour='red', se=FALSE) +
    ylim(0, 1) +
    ylab("cov_95") +
    xlab("log.w") +
    labs(title="95% coverage against Log w, task 4")
  
  multiplot(p_logw68.list[[p]], p_logw95.list[[p]])
  
}