# -----------------------------------------------------------------------------
### Load libraries and source code
#
library(drc)
library(testthat)

# ---------------------------------------------------------------------------
### Do simulation and count the number of times that the Hill bounds do not
### cover the true values of the IC50 and slope paramters.
#
conf.level <- 0.999

n.simul <- 1000
n.obs.dose <- 10
n.settings <- 6  # No. of simulation settings

stdev <- 5

theta.matr <- matrix(c(100, 0.01, -1.25, 0,
                       100, 0.01, -0.75, 0,
                       100, 0.1, -1.25, 0,
                       100, 0.1, -0.75, 0,
                       100, 1, -1.25, 0,
                       100, 1, -0.75, 0),
                     nrow = n.settings,
                     ncol = 4,
                     byrow = TRUE)
levels.dose <- 10^seq(from = -5, to = 0, by = 1)
n.obs <- n.doses*n.obs.dose

result.matr <- matrix(nrow = n.simul, ncol = n.settings)

for(i in 1:n.simul) {

  for(j in 1:n.settings) {

    theta <- theta.matr[j, ]

    x <- rep(levels.dose, each = n.obs.dose)
    y <- MeanResponse(x, theta) + rnorm(length(x), mean = 0, sd = stdev)

    scale.inc <- 0.001
    y.range <- range(y)
    len.y.range <- scale.inc * diff(y.range)

    y.max <- max(y) + len.y.range
    y.min <- min(y) - len.y.range

    y.logit <- log10((y - y.min)/(y.max - y))  # Logit transformed responses
    x.log <- log10(x)  # Log transformed doses

    data.hill <- data.frame(y.logit = y.logit,
                            x.log = x.log)
    data.hill <- subset(data.hill, subset = x.log != -Inf)

    lm.hill <- lm(y.logit ~ x.log, data = data.hill)

    beta <- c(-theta[3]*log10(theta[2]), theta[3])
    beta.bounds <- confint(lm.hill, level = conf.level)

    result.matr[i, j] <- (beta[1] >= beta.bounds[1, 1])&&(beta[1] <= beta.bounds[1, 2])&&
                         (beta[2] >= beta.bounds[2, 1])&&(beta[1] <= beta.bounds[2, 2])
  }
}

colSums(result.matr)
