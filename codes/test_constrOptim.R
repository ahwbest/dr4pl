# ------------------------------------------------------------------------------------
### Load libraries
#
library(drc)

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

FittedValue <- function(theta, x.pred) {

  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]
  
  f <- theta.1 + (theta.4 - theta.1)/(1 + (x.pred/theta.2)^theta.3)
  
  return(f)
}

ErrorFunctionAsTheta1 <- function(theta.1.vec) {
  x.vec <- df.drm$Dose
  y.vec <- df.drm$Response
  
  result <- c()
  
  for(i in 1:length(theta.1.vec)) {
    theta.1 <- theta.1.vec[i]
    
    f.vec <- theta.1 + (theta.4 - theta.1)/(1 + (x.vec/theta.2)^theta.3)
    
    result[i] <- sum((y.vec - f.vec)^2)
  }
  
  return(result)
}

ErrorFunction <- function(theta.vec, df.drm) {
  x.vec <- df.drm$Dose
  y.vec <- df.drm$Response
  
  theta.1 <- theta.vec[1]
  theta.2 <- theta.vec[2]
  theta.3 <- theta.vec[3]
  theta.4 <- theta.vec[4]
  
  f.vec <- theta.1 + (theta.4 - theta.1)/(1 + (x.vec/theta.2)^theta.3)
  
  return(sum((y.vec - f.vec)^2))
}

ErrorFunctionL1 <- function(theta.vec, df.drm) {
  x.vec <- df.drm$Dose
  y.vec <- df.drm$Response
  
  theta.1 <- theta.vec[1]
  theta.2 <- theta.vec[2]
  theta.3 <- theta.vec[3]
  theta.4 <- theta.vec[4]
  
  f.vec <- theta.1 + (theta.4 - theta.1)/(1 + (x.vec/theta.2)^theta.3)
  
  return(sum(abs(y.vec - f.vec)))
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

GradientFunction <- function(theta.vec, df.drm) {
  x.vec <- df.drm$Dose
  y.vec <- df.drm$Response
  
  theta.1 <- theta.vec[1]
  theta.2 <- theta.vec[2]
  theta.3 <- theta.vec[3]
  theta.4 <- theta.vec[4]
  
  eta.vec <- (x.vec/theta.2)^theta.3
  f.vec <- theta.1 + (theta.4 - theta.1)/(1 + eta.vec)
  
  deriv.f.theta.1 <- 1 - 1/(1 + eta.vec)
  deriv.f.theta.2 <- (theta.4 - theta.1)*theta.3/theta.2*eta.vec/(1 + eta.vec)^2
  deriv.f.theta.3 <- -(theta.4 - theta.1)/theta.3*log(eta.vec)*eta.vec/(1 + eta.vec)^2
  deriv.f.theta.4 <- 1/(1 + eta.vec)
  
  return(-2*(y.vec - f.vec)%*%cbind(deriv.f.theta.1, deriv.f.theta.2, deriv.f.theta.3, deriv.f.theta.4))
}

Jacobian <- function(theta.vec, df.drm) {
  x.vec <- df.drm$Dose
  y.vec <- df.drm$Response
  
  theta.1 <- theta.vec[1]
  theta.2 <- theta.vec[2]
  theta.3 <- theta.vec[3]
  theta.4 <- theta.vec[4]
  
  eta.vec <- (x.vec/theta.2)^theta.3
  f.vec <- theta.1 + (theta.4 - theta.1)/(1 + eta.vec)
  
  deriv.f.theta.1 <- 1 - 1/(1 + eta.vec)
  deriv.f.theta.2 <- -(theta.4 - theta.1)*theta.3/theta.2*eta.vec/(1 + eta.vec)^2
  deriv.f.theta.3 <- -(theta.4 - theta.1)/theta.3*log(eta.vec)*eta.vec/(1 + eta.vec)^2
  deriv.f.theta.4 <- 1/(1 + eta.vec)
  
  return(cbind(deriv.f.theta.1, deriv.f.theta.2, deriv.f.theta.3, deriv.f.theta.4))
}

# ------------------------------------------------------------------------------------
### Test the data that works
#

### Load libraries
library(drc)

### Load data
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

df.drm <- read.csv("drc_simulation.csv")

x.vec <- df.drm$Dose
y.vec <- df.drm$Response
n <- length(x.vec)

theta.1 <- max(y.vec)
theta.2 <- median(unique(x.vec))
theta.3 <- -1
theta.4 <- 0.0001

theta.init.vec <- c(theta.1, theta.2, theta.3, theta.4)

constraint.matr <- rbind(c(0, 0, -1, 0), c(0, 0, 0, 1))

drm.cO <- constrOptim(theta = theta.init.vec,
                      f = ErrorFunction, 
                      ui = constraint.matr,
                      ci = c(0, 0),
                      df.drm = df.drm,
                      method = "Nelder-Mead")

theta.fitted.vec <- drm.cO$par
theta.1 <- theta.fitted.vec[1]
theta.2 <- theta.fitted.vec[2]
theta.3 <- theta.fitted.vec[3]
theta.4 <- theta.fitted.vec[4]

f.vec <- theta.1 + (theta.4 - theta.1)/(1 + (x.vec/theta.2)^theta.3)

jacobian <- Jacobian(theta.fitted.vec, df.drm)

C.hat.inv <- solve(t(jacobian)%*%jacobian)
s <- sqrt(sum((y.vec - f.vec)^2)/(n - 4))

q.t <- qt(0.975, df = n - 4)
ci <- cbind(theta.fitted.vec - q.t*s*sqrt(diag(C.hat.inv)), theta.fitted.vec + q.t*s*sqrt(diag(C.hat.inv)))
ci <- round(ci, 3)

drm.drc <- drm(Response ~ Dose,
             data = df.drm,
             fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50")))

plot(x = log10(x.vec),
     y = y.vec,
     main = "Scatter plot",
     xlab = "",
     ylab = "Response",
     lwd = 2,
     pch = 16,
     type = "p",
     bty = "L",
     xaxt = "n")

theta.vec = theta.fitted.vec
curve(DoseResponseCurve, type = "l", add = TRUE, lwd = 2)
mtext(text = "Dose", side = 1, line = 4)

x.levels <- sort(unique(x.vec))
at.x <- log10(x.levels)
axis(1, at = at.x, labels = 10^at.x,las=2)

# ------------------------------------------------------------------------------------
### Test the data that fails
#



### Set the output file
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Presentation")
setEPS()
postscript("constrOptim_gof_3.eps", width = 4, height = 4)

### Load data
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

par(mar = c(5, 4, 2, 1) + 0.1)

df.drm <- read.csv("drc_error_3.csv")

x <- df.drm$Dose
y <- df.drm$Response

theta.1 <- max(y)
theta.2 <- median(unique(x))
theta.3 <- -1
theta.4 <- min(y) + 0.0001

theta.init <- c(theta.1, theta.2, theta.3, theta.4)

constraint.matr <- rbind(c(0, 0, -1, 0), c(0, 0, 0, 1))

drm.cO <- constrOptim(theta = theta.init,
                      f = ErrorFunction, 
                      ui = constraint.matr,
                      ci = c(0, 0),
                      df.drm = df.drm,
                      method = "Nelder-Mead")

theta.fitted <- drm.cO$par
gof.pval <- GoodnessOfFit(theta.fitted, df.drm)

plot(x = log10(x),
     y = y,
     main = "Scatter plot (L2-error)",
     xlab = "",
     ylab = "Response",
     lwd = 2,
     pch = 16,
     type = "p",
     bty = "L",
     xaxt = "n")

theta.l2 = theta.fitted
theta = theta.l2

curve(DoseResponseCurve, type = "l", add = TRUE, lwd = 2)
mtext(text = "Dose", side = 1, line = 4)

legend(x = "topright", bty = "n", legend = paste("gof p-val =", round(gof.pval,3)))

x.levels <- sort(unique(x))
at.x <- log10(x.levels)
axis(1, at = at.x, labels = 10^at.x,las=2)

dev.off()