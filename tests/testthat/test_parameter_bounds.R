# -----------------------------------------------------------------------------
### Load libraries
#
library(drc)
library(testthat)

# ---------------------------------------------------------------------------
### Do simulation and count the number of times that the Hill bounds do not
### cover the true values of the IC50 and slope paramters.
#
conf.level <- 0.9999
stdev <- 5

n.simul <- 1000
n.obs.dose <- 5
n.settings <- 6  # No. of simulation settings

theta.matr <- matrix(c(100, 0.01, -1.25, 0,
                       100, 0.01, -0.75, 0,
                       100, 0.1, -1.25, 0,
                       100, 0.1, -0.75, 0,
                       100, 1, -1.25, 0,
                       100, 1, -0.75, 0),
                     nrow = n.settings,
                     ncol = 4,
                     byrow = TRUE)
levels.dose <- 10^seq(from = -4, to = 1, by = 1)
n.doses <- length(levels.dose)
n.obs <- n.doses*n.obs.dose

result.matr <- matrix(nrow = n.simul, ncol = n.settings)

for(i in 1:n.simul) {

  for(j in 1:n.settings) {

    theta <- theta.matr[j, ]

    x <- rep(levels.dose, each = n.obs.dose)
    y <- MeanResponse(x, theta) + rnorm(length(x), mean = 0, sd = stdev)

    theta.init <- FindInitialParms(x, y, "Mead", NULL)
    
    theta.init.1 <- theta.init[1]
    theta.init.4 <- theta.init[4]
    
    data.hill <- data.frame(x = x, y = y)
    data.hill <- subset(data.hill, subset = (data.hill$y<theta.init.1)&
                                            (data.hill$y>theta.init.4)&
                                            (data.hill$x>0))
    
    data.hill$y.logit <- log((data.hill$y - theta.init.4)/(theta.init.1 - data.hill$y))  # Logit transformed responses
    data.hill$x.log <- log(data.hill$x)  # Log transformed doses
    
    lm.hill <- lm(y.logit ~ x.log, data = data.hill)
    
    beta.bounds <- confint(lm.hill, level = 0.999)
    beta.2 <- -theta[3]*log(theta[2])
    beta.3 <- theta[3]
    
    result.matr[i, j] <- (beta.2 >= beta.bounds[1, 1])&&(beta.2 <= beta.bounds[1, 2])&&
                         (beta.3 >= beta.bounds[2, 1])&&(beta.3 <= beta.bounds[2, 2])
  }
}

colSums(result.matr)
