
library(drc)

# --------------------------------------------------------------------------------
### User-defined functions
#
ErrorForTheta2and3 <- function(theta2.vec, theta3.vec) {
  # Need x.vec and y.vec in advance
  n.row <- length(theta2.vec)
  n.col <- length(theta3.vec)
  error.matr <- matrix(data = 0, nrow = n.row, ncol = n.col)
  
  for(i in 1:n.row) {
    for(j in 1:n.col) {
      theta2 <- theta2.vec[i]
      theta3 <- theta3.vec[j]
      
      f.vec <- theta1 + (theta4 - theta1)/(1 + 10^(theta3*(log10(x.vec) - theta2)))
      error.matr[i, j] <- sum((y.vec - f.vec)^2)
    }
  }
  
  return(error.matr)
}

# --------------------------------------------------------------------------------
### Analysis
#
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

df.drc <- read.csv("drc_error_2.csv")
df.drc <- df.drc[df.drc$Dose != 0, ]

x.orig.vec <- df.drc$Dose
y.orig.vec <- df.drc$Response

### Scatter plot and contour plot
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Presentation")
setEPS()
postscript("error_function_plot_failure_1_2.eps", width = 8, height = 4)

par(mfrow = c(1, 2), mar = c(5, 4, 3, 2) + 0.1)

# Scatter plot
plot(x = log10(x.orig.vec),
     y = y.orig.vec,
     main = "Scatter plot",
     xlab = "",
     ylab = "Response",
     lwd = 2,
     pch = 16,
     type = "p",
     bty = "L",
     xaxt = "n")

mtext(text = "Dose", side = 1, line = 4)

x.levels <- sort(unique(x.orig.vec))
at.x <- log10(x.levels)
axis(1, at = at.x, labels = 10^at.x, las = 2)

# Contour plot
x.vec <- x.orig.vec
y.vec <- y.orig.vec

theta1 <- max(y.vec)
theta4 <- min(y.vec)

theta2.vec <- seq(from = log10(min(x.vec)), to = log10(max(x.vec)), length = 100)
theta3.vec <- seq(from = -1, to = 0, length = 100)
error.matr <- ErrorForTheta2and3(theta2.vec, theta3.vec)

contour(x = theta2.vec,
        y = theta3.vec,
        z = error.matr,
        main = "Contour plot",
        xlab = "",
        ylab = "Theta3",
        nlevels = 10,
        labcex = 1,
        xaxt = "n",
        bty = "l")
at.x=seq(from = min(theta2.vec), to = max(theta2.vec), by = 1)
#at.x=seq(from = min(floor(log10(x.vec))), to = max(floor(log10(x.vec))), by = 1)
axis(1, at = at.x, labels = 10^at.x, las=2)
mtext(text = "Theta2", side = 1, line = 4)

dev.off()