# -----------------------------------------------------------------------------
### Load libraries
#
library(drc)

# ----------------------------------------------------------------------------
### Count the numbers of times that the confidence intervals contain the true
### parameters
#
conf.level <- 0.9999
stdev <- 5

n.simul <- 1000
n.obs.dose <- 5
n.settings <- 9  # No. of simulation settings

theta.mat <- matrix(c(100, 0.01, -0.5, 0,
                      100, 0.01, -1, 0,
                      100, 0.01, -1.5, 0,
                      100, 0.1, -0.5, 0,
                      100, 0.1, -1, 0,
                      100, 0.1, -1.5, 0,
                      100, 1, -0.5, 0,
                      100, 1, -1, 0,
                      100, 1, -1.5, 0),
                    nrow = n.settings,
                    ncol = 4,
                    byrow = TRUE)
levels.dose <- 10^seq(from = -4, to = 0, by = 1)
n.doses <- length(levels.dose)
n.obs <- n.doses*n.obs.dose

result.mat <- matrix(nrow = n.simul, ncol = n.settings)

for(i in 1:n.simul) {
  
  for(j in 1:n.settings) {
    
    theta <- theta.mat[j, ]
    
    x <- rep(levels.dose, each = n.obs.dose)
    y <- MeanResponse(x, theta) + rnorm(length(x), mean = 0, sd = stdev)
    
    ## Obtain initial parameter estimates.
    theta.init <- FindInitialParms(x, y,
                                   method.init = "Mead")
    
    retheta.init <- ParmToLog(theta.init)  # Reparameterized parameters
    Hill.bounds <- FindHillBounds(x, y, retheta.init, use.Hessian = FALSE)
    
    result.mat[i, j] <- (log10(theta[2]) >= Hill.bounds$LogTheta2[1])&&
                        (log10(theta[2]) <= Hill.bounds$LogTheta2[2])&&
                        (theta[3] >= Hill.bounds$Theta3[1])&&
                        (theta[3] <= Hill.bounds$Theta3[2])
  }
}

results <- rep(0, ncol(result.mat))

for(j in 1:ncol(result.mat)) {
  
  results[j] <- mean(result.mat[!is.na(result.mat[, j]), j])
  
}

cat(1 - results)

# ----------------------------------------------------------------------------
### Count the numbers of times that the confidence intervals are attained
### during optimization processes.
#
conf.level <- 0.9999
stdev <- 5

n.simul <- 1000
n.obs.dose <- 5
n.settings <- 9  # No. of simulation settings

theta.mat <- matrix(c(100, 0.01, -0.5, 0,
                      100, 0.01, -1, 0,
                      100, 0.01, -1.5, 0,
                      100, 0.1, -0.5, 0,
                      100, 0.1, -1, 0,
                      100, 0.1, -1.5, 0,
                      100, 1, -0.5, 0,
                      100, 1, -1, 0,
                      100, 1, -1.5, 0),
                    nrow = n.settings,
                    ncol = 4,
                    byrow = TRUE)
levels.dose <- 10^seq(from = -4, to = 0, by = 1)
n.doses <- length(levels.dose)
n.obs <- n.doses*n.obs.dose

result.mat <- matrix(nrow = n.simul, ncol = n.settings)

for(i in 1:n.simul) {
  
  for(j in 1:n.settings) {
    
    theta <- theta.mat[j, ]
    
    x <- rep(levels.dose, each = n.obs.dose)
    y <- MeanResponse(x, theta) + rnorm(length(x), mean = 0, sd = stdev)
    
    data.simul <- data.frame(Dose = x, Response = y)
    
    dr4pl.simul <- dr4pl(Response ~ Dose,
                         data = data.simul)

    result.mat[i, j] <- dr4pl.simul$convergence
    
  }
  
}

results <- rep(0, ncol(result.mat))

for(j in 1:ncol(result.mat)) {
  
  results[j] <- mean(result.mat[!is.na(result.mat[, j]), j])
  
}

cat(results)