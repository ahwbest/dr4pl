### Load libraries
library(drc)

### User-defined functions
DoseResponseCurve <- function(log10.x.vec) {
  x.vec <- 10^log10.x.vec
  
  f.vec <- theta.1 + (theta.4 - theta.1)/(1 + (x.vec/theta.2)^theta.3)
  
  return(f.vec)
}

DRC5PL <- function(log10.x.vec) {
  x.vec <- 10^log10.x.vec
  
  f.vec <- theta.1 + (theta.4 - theta.1)/(1 + (x.vec/theta.2)^theta.3)^theta.5
  
  return(f.vec)
}

ErrorForTheta2 <- function(theta2.vec) {
  # Need x.vec and y.vec in advance
  theta1 <- max(y.vec)
  theta3 <- -1
  theta4 <- min(y.vec)
  
  error.vec <- rep(0, length(theta2.vec))
  for(i in 1:length(theta2.vec)) {
    theta2 <- theta2.vec[i]
    
    f.vec <- theta1 + (theta4 - theta1)/(1 + 10^(theta3*(log10(x.vec) - log10(theta2))))
    
    error.vec[i] <- sum((y.vec - f.vec)^2)
  }
  
  return(error.vec)
}

# -----------------------------------------------------------------------------------
### Simulation
#
theta.1 <- 100
theta.2 <- 1.35
theta.3 <- -1
theta.4 <- 0
theta.5 <- 1.2

x.levels <- c(0.00135, 0.0135, 0.135, 1.35, 13.5, 135, 1350)
x.vec <- rep(x.levels, each = 5)
f.vec <- theta.1 + (theta.4 - theta.1)/(1 + 10^(theta.3*(log10(x.vec) - log10(theta.2))))
y.vec <- f.vec + rnorm(n = length(x.vec), mean = 0, sd = sqrt(400))

# --------------------------------------------------------------------------------------
### Plot a dose response curve
#
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Presentation")
setEPS()
postscript("4PL_5PL.eps", width = 6, height = 4)

curve(DoseResponseCurve, type = "l", lwd = 2,
      xlab = "",
      ylab = "Response",
      bty = "l",
      xaxt = "n",
      xlim = c(log10(0.00135), log10(1350)))

curve(DRC5PL, type = "l", lwd = 2, add = TRUE, col = "red")

mtext(text = "Dose", side = 1, line = 4)

at.x=log10(x.levels)
axis(1,at=at.x,labels=10^at.x,las=2)

legend(x = "topright", legend = c("4PL", "5PL"), 
       col = c("black", "red"),
       lwd = 2,
       bty = "n")

dev.off()

# --------------------------------------------------------------------------------------
### Plot simulated data
#
setwd("C:\\Users\\Hyowon\\Presentation\\Dose_response_modelling")
setEPS()
postscript("simulated_data.eps", width = 6, height = 4)

plot(drm.plot,
     xlab = "Dose",
     ylab = "Response",
     lwd = 2,
     pch = 16,
     type = "all",
     bty = "L")
curve(DoseResponseCurve, type = "l", add = TRUE, lwd = 2,
      col = "red")
legend(x = "topright", legend = c("True mean", "Fitted mean"), 
       col = c("black", "red"),
       lwd = 2,
       bty = "n")

dev.off()

# -------------------------------------------------------------------------------
### Save successful cases as csv files
#
df.simul <- as.data.frame(cbind(x.vec.orig, y.vec.orig))
colnames(df.simul) <- c("Dose", "Response")

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")
write.csv(df.simul, "drc_simulation1.csv", row.names = FALSE)

### Load x.vec.orig and y.vec.orig here

# -------------------------------------------------------------------------------
### Remove right most points sequentially
#
length.x.vec <- c(35, 30, 25, 20, 15, 10)

setwd("C:\\Users\\Hyowon\\Presentation\\Dose_response_modelling")
setEPS()
postscript("removing_points.eps", width = 6, height = 4)

par(mar = c(6,4,2,2)+0.1, mfrow = c(2, 3))
for(length.x in length.x.vec) {
  x.vec <- x.vec.orig[1:length.x]
  y.vec <- y.vec.orig[1:length.x]
  
  df.drm <- as.data.frame(cbind(x.vec, y.vec))
  colnames(df.drm) <- c("Dose", "Response")
  
  drm.simul <- tryCatch(drm(Response ~ Dose,
                            data = df.drm,
                            fct = LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "IC50"))),
                        error = function(e) { NULL })
  
  if(class(drm.simul) == "drc") {
    plot(drm.simul,
         main = paste("n =", length.x),
         xlab = "",
         ylab = "Response",
         xlim = c(min(x.vec.orig), max(x.vec.orig)),
         ylim = c(min(y.vec.orig), max(y.vec.orig)),
         lwd = 2,
         legend = FALSE,
         pch = 16,
         type = "all",
         bty = "l")
    
    mtext(text = "Dose", side = 1, line = 4, cex = 0.7)
  } else {
    plot(x = log10(x.vec),
         y = y.vec,
         main = paste("n =", length.x),
         xlab="",
         ylab="Response",
         xlim = c(min(log10(x.vec.orig)), max(log10(x.vec.orig))),
         ylim = c(min(y.vec.orig), max(y.vec.orig)),
         lwd=2,
         pch=16,
         type="p",
         bty="L",
         xaxt = "n")

    mtext(text = "Dose", side = 1, line = 4, cex = 0.7)
    
    at.x=log10(x.levels)
    axis(1,at=at.x,labels=10^at.x,las=2)
  }
}

dev.off()