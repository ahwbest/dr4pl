# -------------------------------------------------------------------------------------
### User-defined functions
#
DoseResponseCurve <- function(log10.x.vec) {
  x.vec <- 10^log10.x.vec
  
  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]
  
  f.vec <- theta.1 + (theta.4 - theta.1)/(1 + (x.vec/theta.2)^theta.3)
  
  return(f.vec)
}

DRC5PL <- function(log10.x.vec) {
  x.vec <- 10^log10.x.vec
  
  f.vec <- theta.1 + (theta.4 - theta.1)/(1 + (x.vec/theta.2)^theta.3)^theta.5
  
  return(f.vec)
}

ErrorFunction <- function(theta, df.drm) {
  x <- df.drm$Dose
  y <- df.drm$Response
  
  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]
  
  f <- theta.1 + (theta.4 - theta.1)/(1 + (x/theta.2)^theta.3)
  
  return(sum((y - f)^2))
}

FittedValue <- function(theta, x.pred) {
  
  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]
  
  f <- theta.1 + (theta.4 - theta.1)/(1 + (x.pred/theta.2)^theta.3)
  
  return(f)
}

GoodnessOfFit <- function(theta.fitted, df.drm) {
  
  x <- df.drm$Dose
  y <- df.drm$Response
  
  x.level <- sort(unique(x))
  y.bar <- tapply(X = y, INDEX = x, FUN = mean)
  mu.hat <- FittedValue(theta.fitted, x.level)
  
  J.i <- table(x)
  n <- length(unique(x))
  N <- length(x)
  p <- length(theta.fitted)
  
  y.bar.vec <- c()
  x.idx <- as.numeric(names(y.bar))
  
  if(all(x.idx == x.level) == FALSE) {
    stop("The order of dose levels has been changed while be executed.")
  }
  
  for(i in 1:nrow(df.drm)) {
    y.bar.vec[i] <- y.bar[which(x.idx == x[i])]
  }
  
  gof.numer <- J.i%*%(y.bar - mu.hat)^2/(n - p)
  gof.denom <- sum((y - y.bar.vec)^2)/(N - n)
  
  gof.pval <- pf(gof.numer/gof.denom, df1 = n - p, df2 = N - n, lower.tail = FALSE)
  
  return(gof.pval)
}

# ---------------------------------------------------------------------------
### Simulation
#

# Set these parameters before running the following codes
n.rep.vec <- c(2, 3, 10)
n.simul <- 10000
sigma <- 20

x.level <- c(0.00135, 0.0135, 0.135, 1.35, 13.5)

theta.1 <- 80
theta.2 <- 0.135
theta.3 <- -1
theta.4 <- 20

result <- matrix(0, nrow = length(n.rep.vec), 3)

for(r in 1:length(n.rep.vec)) {
  n.rep <- n.rep.vec[r]
  pval <- rep(0, n.simul)
  
  for(i in 1:n.simul) {
    x <- rep(x.level, each = n.rep)
    f <- theta.1 + (theta.4 - theta.1)/(1 + (x/theta.2)^theta.3)
    y <- f + rnorm(n = length(x), mean = 0, sd = sigma)
    df.drm <- data.frame(Dose = x, Response = y)
    
    theta.init <- c(theta.1, theta.2, theta.3, theta.4)
    constr.matr <- c(0, 0, -1, 0)
    
    drm.cO <- constrOptim(theta = theta.init,
                          f = ErrorFunction, 
                          ui = constr.matr,
                          ci = 0,
                          df.drm = df.drm,
                          method = "Nelder-Mead")
    
    theta.fitted <- drm.cO$par
    pval[i] <- GoodnessOfFit(theta.fitted, df.drm)
  }
  
  result[r, ] <- c(sum(pval < 0.1)/n.simul, sum(pval < 0.05)/n.simul, sum(pval < 0.025)/n.simul)
}

print(round(result, 3))

theta = c(theta.1, theta.2, theta.3, theta.4)

### Set the output file
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Presentation")
setEPS()
postscript("gof_size_2.eps", width = 6, height = 4)

curve(DoseResponseCurve,
      main = "Scatter plot (L2-error)",
      xlab = "",
      ylab = "Response",
      xlim = c(min(log10(x.level)), max(log10(x.level))),
      ylim = c(0, 100),
      type = "l",
      bty = "l",
      xaxt = "n",
      lwd = 2)
mtext(text = "Dose", side = 1, line = 4)

x.levels <- sort(unique(x))
at.x <- log10(x.levels)
axis(1, at = at.x, labels = 10^at.x,las=2)

dev.off()

# ---------------------------------------------------------------------------
### Simulation for the departure from a null model
#

# Set these parameters before running the following codes
n.rep.vec <- c(300)
n.simul <- 10000
sigma <- 10

x.level <- c(0.00135, 0.0135, 0.135, 1.35, 13.5)

theta.1 <- 100
theta.2 <- 0.135
theta.3 <- -1
theta.4 <- 0
theta.5 <- 1.2

result <- matrix(0, nrow = length(n.rep.vec), 3)

for(r in 1:length(n.rep.vec)) {
  n.rep <- n.rep.vec[r]
  pval <- rep(0, n.simul)
  
  for(i in 1:n.simul) {
    x <- rep(x.level, each = n.rep)
    f <- theta.1 + (theta.4 - theta.1)/(1 + (x/theta.2)^theta.3)^theta.5
    y <- f + rnorm(n = length(x), mean = 0, sd = sigma)
    df.drm <- data.frame(Dose = x, Response = y)
    
    theta.init <- c(theta.1, theta.2, theta.3, theta.4)
    constr.matr <- c(0, 0, -1, 0)
    
    drm.cO <- constrOptim(theta = theta.init,
                          f = ErrorFunction, 
                          ui = constr.matr,
                          ci = 0,
                          df.drm = df.drm,
                          method = "Nelder-Mead")
    
    theta.fitted <- drm.cO$par
    pval[i] <- GoodnessOfFit(theta.fitted, df.drm)
  }
  
  result[r, ] <- c(sum(pval < 0.1)/n.simul, sum(pval < 0.05)/n.simul, sum(pval < 0.025)/n.simul)
}

print(round(result, 3))

theta = c(theta.1, theta.2, theta.3, theta.4)

### Set the output file
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Presentation")
setEPS()
postscript("gof_size_2.eps", width = 6, height = 4)

curve(DoseResponseCurve,
      main = "Scatter plot (L2-error)",
      xlab = "",
      ylab = "Response",
      xlim = c(min(log10(x.level)), max(log10(x.level))),
      ylim = c(0, 100),
      type = "l",
      bty = "l",
      xaxt = "n",
      lwd = 2)
mtext(text = "Dose", side = 1, line = 4)

x.levels <- sort(unique(x))
at.x <- log10(x.levels)
axis(1, at = at.x, labels = 10^at.x,las=2)

dev.off()