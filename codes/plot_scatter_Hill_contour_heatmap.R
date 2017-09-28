
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

df.drc <- read.csv("drc_simulation.csv")
df.drc <- df.drc[df.drc$Dose != 0,]

x.vec.orig <- df.drc$Dose
y.vec.orig <- df.drc$Response


setwd("C:\\Users\\Hyowon\\Presentation\\Dose_response_modelling")
setEPS()
postscript("drc_failure_simulation_1.eps", width = 4, height = 4)

plot(x = log10(df.drc$Dose),
     y = df.drc$Response,
     main = "Scatter plot",
     xlab = "",
     ylab = "Response",
     lwd = 2,
     pch = 16,
     type = "p",
     bty = "L",
     xaxt = "n")

mtext(text = "Dose", side = 1, line = 4)

x.levels <- sort(unique(x.vec.orig))
at.x=log10(x.levels)
axis(1,at=at.x,labels=10^at.x,las=2)

dev.off()

### Hill plot
setEPS()
postscript("drc_failure_simulation_2.eps", width = 4, height = 4)

x.vec <- x.vec.orig
y.vec <- y.vec.orig

y.transformed.vec <- log10((y.vec-min(y.vec))/(max(y.vec)-y.vec))
idx.non.inf <- which(y.transformed.vec == Inf | y.transformed.vec == -Inf)
y.transformed.vec <- y.transformed.vec[-idx.non.inf]
x.vec <- x.vec.orig[-idx.non.inf]
lm.test <- lm(y.transformed.vec ~ log10(x.vec))
R2 <- summary(lm.test)$r.squared

plot(x = log10(x.vec), y = y.transformed.vec, 
     main = "Hill plot",
     xlab = "",
     ylab = "Transformed response",
     bty = "l",
     type = "p",
     pch = 16,
     xaxt = "n")

legend(x = "topright", legend = paste("R2 =", round(R2,3)), bty = "n")

mtext(text = "Dose", side = 1, line = 4)

x.levels <- sort(unique(x.vec.orig))
at.x=log10(x.levels)
axis(1,at=at.x,labels=10^at.x,las=2)

dev.off()

### Contour plot
setEPS()
postscript("drc_failure_simulation_3.eps", width = 4, height = 4)

x.vec <- x.vec.orig
y.vec <- y.vec.orig

theta1 <- max(y.vec)
theta4 <- min(y.vec)

theta2.vec <- seq(from = log10(min(x.vec)), to = log10(max(x.vec)), length = 100)
theta3.vec <- seq(from = -1, to = 0, length = 100)
error.matr <- ErrorForTheta2and3(theta2.vec, theta3.vec)

contour(x = theta2.vec,
        y = theta3.vec,
        z = error.matr,
        main = "Contour plot",
        xlab = "Theta2",
        ylab = "Theta3",
        nlevels = 10,
        labcex = 1,
        xaxt = "n")
at.x=seq(from = min(floor(log10(x.vec))), to = max(floor(log10(x.vec))), by = 1)
axis(1, at = at.x, labels = 10^at.x, las=2)

dev.off()

### Heatmap
setEPS()
postscript("drc_failure_simulation_4.eps", width = 4, height = 4)

x.vec <- x.vec.orig
y.vec <- y.vec.orig

theta1 <- max(y.vec)
theta4 <- min(y.vec)

theta2.vec <- seq(from = log10(min(x.vec)), to = log10(max(x.vec)), length = 100)
theta3.vec <- seq(from = -1, to = 0, length = 100)
error.matr <- ErrorForTheta2and3(theta2.vec, theta3.vec)

heatmap(x = t(log10(error.matr)), Rowv = NA, Colv = NA,
        labRow = "",
        labCol = "",
        scale = "none")

dev.off()
