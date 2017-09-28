# ------------------------------------------------------------------------
### Load libraries
#
library(drc)

# --------------------------------------------------------------------------------
### User-defined functions
#
ErrorForBeta2and3 <- function(beta2.vec, beta3.vec) {

  # Need x and y in advance
  n.row <- length(beta2.vec)
  n.col <- length(beta3.vec)
  error.matr <- matrix(data = 0, nrow = n.row, ncol = n.col)

  for(i in 1:n.row) {

    for(j in 1:n.col) {

      beta2 <- beta2.vec[i]
      beta3 <- beta3.vec[j]

      f <- theta1 + (theta4 - theta1)/(1 + 10^(theta3*(log10(x) - theta2)))
      error.matr[i, j] <- sum((y - f)^2)/length(y)
    }
  }

  return(error.matr)
}

ErrorForTheta2and3 <- function(theta2.vec, theta3.vec) {

  # Need x and y in advance
  n.row <- length(theta2.vec)
  n.col <- length(theta3.vec)
  error.matr <- matrix(data = 0, nrow = n.row, ncol = n.col)

  for(i in 1:n.row) {

    for(j in 1:n.col) {

      theta2 <- theta2.vec[i]
      theta3 <- theta3.vec[j]

      f <- theta1 + (theta4 - theta1)/(1 + 10^(theta3*(log10(x) - theta2)))
      error.matr[i, j] <- sum((y - f)^2)/length(y)
    }
  }

  return(error.matr)
}

# ---------------------------------------------------------------------------
### Load data
#
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

data.error.list <- vector("list", 4)

data.error.1 <- read.csv("drc_error_1.csv")
data.error.2 <- read.csv("drc_error_2.csv")
data.error.3 <- read.csv("drc_error_3.csv")
data.error.4 <- read.csv("drc_error_4.csv")

for(i in 1:4) {

  file.name <- paste("drc_error_", i, ".csv", sep = "")
  data.error.list[[i]] <- read.csv(file.name)
}

# ---------------------------------------------------------------------------
### Draw scatter plots
#
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")
jpeg(filename = "drc_error_cases.jpg",
     width = 8, height = 8, units = "in", res = 600)

par(mfrow = c(2, 2), mar = c(6, 4, 3, 1) + 0.1)

for(i in 1:4) {

  data.error <- data.error.list[[i]]
  doses <- data.error$Dose
  responses <- data.error$Response

  plot(x = log10(doses), y = responses,
       pch = 16,
       main = paste("Error Case ", i),
       xlab = "",
       ylab = "Response",
       bty = "l",
       xaxt = "n")
  axis(side = 1, at = log10(unique(doses)),
       labels = unique(doses),
       las = 2)
  mtext(side = 1, text = "Dose", line = 4)
}

dev.off()

### Draw Hill plots

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")
jpeg(filename = "Hill_plot_error_cases.jpg",
     width = 8, height = 8, units = "in", res = 600)

par(mfrow = c(2, 2), mar = c(5, 4, 3, 1) + 0.1)

x = data.error.1$Dose
y = data.error.1$Response

scale.inc <- 0.001
y.range <- range(y)
len.y.range <- scale.inc * diff(y.range)

y.max <- max(y) + len.y.range
y.min <- min(y) - len.y.range

y.logit <- log10((y - y.min)/(y.max - y))  # Logit transformed responses
x.log <- log10(x)  # Log transformed doses

plot(x = x.log, y = y.logit, pch = 16,
     main = "Error Case 1",
     xlab = expression("log"[10]*"(dose)"),
     ylab = "logit(Response)",
     bty = "l",
     xaxt = "n")

axis(side = 1, at = log10(unique(data.error.1$Dose)),
     labels = unique(data.error.1$Dose))

x <- data.error.2$Dose
y <- data.error.2$Response

scale.inc <- 0.001
y.range <- range(y)
len.y.range <- scale.inc * diff(y.range)

y.max <- max(y) + len.y.range
y.min <- min(y) - len.y.range

y.logit <- log10((y - y.min)/(y.max - y))  # Logit transformed responses
x.log <- log10(x)  # Log transformed doses

plot(x = x.log, y = y.logit, pch = 16,
     main = "Error Case 2",
     xlab = expression("log"[10]*"(dose)"),
     ylab = "logit(Response)",
     bty = "l",
     xaxt = "n")

axis(side = 1, at = log10(unique(data.error.2$Dose)),
     labels = unique(data.error.2$Dose))

x <- data.error.3$Dose
y <- data.error.3$Response

scale.inc <- 0.001
y.range <- range(y)
len.y.range <- scale.inc * diff(y.range)

y.max <- max(y) + len.y.range
y.min <- min(y) - len.y.range

y.logit <- log10((y - y.min)/(y.max - y))  # Logit transformed responses
x.log <- log10(x)  # Log transformed doses

plot(x = x.log, y = y.logit, pch = 16,
     main = "Error Case 3",
     xlab = expression("log"[10]*"(dose)"),
     ylab = "logit(Response)",
     bty = "l",
     xaxt = "n")

axis(side = 1, at = log10(unique(data.error.3$Dose)),
     labels = unique(data.error.3$Dose))

x <- data.error.4$Dose
y <- data.error.4$Response

scale.inc <- 0.001
y.range <- range(y)
len.y.range <- scale.inc * diff(y.range)

y.max <- max(y) + len.y.range
y.min <- min(y) - len.y.range

y.logit <- log10((y - y.min)/(y.max - y))  # Logit transformed responses
x.log <- log10(x)  # Log transformed doses

plot(x = x.log, y = y.logit, pch = 16,
     main = "Error Case 4",
     xlab = expression("log"[10]*"(dose)"),
     ylab = "logit(Response)",
     bty = "l",
     xaxt = "n")

axis(side = 1, at = log10(unique(data.error.4$Dose)),
     labels = unique(data.error.4$Dose))

dev.off()

# ----------------------------------------------------------------------------
### Draw contour plots
#
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")
jpeg(filename = "contour_plot_error_cases.jpg",
     width = 8, height = 8, units = "in", res = 600)

par(mfrow = c(2, 2), mar = c(6, 4, 2, 1) + 0.1)

for(i in 1:4) {

  data.error <- data.error.list[[i]]
  doses <- data.error$Dose
  responses <- data.error$Response

  x <- doses[doses != 0]
  y <- responses[doses != 0]

  scale.inc <- 0.001
  y.range <- range(y)
  len.y.range <- scale.inc * diff(y.range)

  theta1 <- max(y) + len.y.range
  theta4 <- min(y) - len.y.range

  beta0.vec <- seq(from = log10(min(x)), to = log10(max(x)), length = 300)
  beta1.vec <- seq(from = -1, to = 0, length = 300)
  error.matr <- ErrorForTheta2and3(beta0.vec, beta1.vec)

  contour(x = beta0.vec,
          y = beta1.vec,
          z = error.matr,
          main = "",
          xlab = "",
          ylab = expression(theta[3]),
          nlevels = 20,
          labcex = 1,
          xaxt = "n",
          bty = "l",
          cex.lab = 1.2)

  at.x=seq(from = floor(log10(min(x))), to = ceiling(log10(max(x))), by = 1)
  axis(1, at = at.x, labels = 10^(at.x), las=2)

  mtext(text = paste("Error Case ", i), side = 3, line = 1, cex = 1.2)
  mtext(text = expression(theta[2]), side = 1, line = 4)
}

dev.off()

### Draw contour plots with reparametrization

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")
jpeg(filename = "reparametrized_contour_plot_error_cases.jpg",
     width = 8, height = 4, units = "in", res = 600)

par(mfrow = c(1, 2), mar = c(5, 4, 1, 1) + 0.1)

x <- data.error.1$Dose[data.error.1$Dose != 0]
y <- data.error.1$Response[data.error.1$Dose != 0]

scale.inc <- 0.001
y.range <- range(y)
len.y.range <- scale.inc * diff(y.range)

theta1 <- max(y) + len.y.range
theta4 <- min(y) - len.y.range

beta0.vec <- seq(from = log10(min(x)), to = log10(max(x)), length = 300)
beta1.vec <- seq(from = -1, to = 0, length = 300)
error.matr <- ErrorForTheta2and3(beta0.vec, beta1.vec)

contour(x = beta0.vec,
        y = beta1.vec,
        z = error.matr,
        main = "Error Case 1",
        xlab = "",
        ylab = expression(theta[3]),
        nlevels = 20,
        labcex = 1,
        xaxt = "n",
        bty = "l")

at.x=seq(from = floor(log10(min(x))), to = ceiling(log10(max(x))), by = 1)
axis(1, at = at.x, labels = 10^(at.x), las=2)
mtext(text = expression(theta[2]), side = 1, line = 4)

x <- data.error.2$Dose[data.error.2$Dose != 0]
y <- data.error.2$Response[data.error.2$Dose != 0]

scale.inc <- 0.001
y.range <- range(y)
len.y.range <- scale.inc * diff(y.range)

theta1 <- max(y) + len.y.range
theta4 <- min(y) - len.y.range

beta0.vec <- seq(from = log10(min(x)), to = log10(max(x)), length = 300)
beta1.vec <- seq(from = -1, to = 0, length = 300)
error.matr <- ErrorForTheta2and3(beta0.vec, beta1.vec)

contour(x = beta0.vec,
        y = beta1.vec,
        z = error.matr,
        main = "Error Case 2",
        xlab = "",
        ylab = expression(theta[3]),
        nlevels = 20,
        labcex = 1,
        xaxt = "n",
        bty = "l")

at.x <- log10(unique(x))
axis(1, at = at.x, labels = 10^(at.x), las=2)
mtext(text = expression(theta[2]), side = 1, line = 4)

dev.off()
