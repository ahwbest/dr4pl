# --------------------------------------------------------------------------
### User-defined functions for M-estimators
#
squared.loss <- function(r) {
  return(r^2)
}

absolute.loss <- function(r) {
  return(abs(r))
}

huber.loss <- function(r) {
  # Values of Huber's loss function evaluated at residuals r.
  #
  # Args:
  #   r: Residuals
  #
  # Returns:
  #   result: Huber's loss function values evaluated at residuals r.
  const <- 1.345  # This value should be chosen in an adaptive fashion.

  ret.val <- r^2  # Vector of results
  outer.term <- 2*const*abs(r) - const^2

  outer.idx <- abs(r) > const

  ret.val[outer.idx] <- outer.term[outer.idx]

  return(ret.val)
}

tukey.biweight.loss <- function(r) {
  # Values of Tukey's biweight loss function evaluated at residuals r.
  #
  # Args:
  #   r: Residuals
  #
  # Returns:
  #   result: Tukey's biweight loss function values evaluated at residuals r.
  const <- 4.685

  ret.val <- (r^6)/(const^4) - 3*(r^4)/(const^2) + 3*r^2

  ret.val[abs(r) > const] <- const^2

  return(ret.val)
}

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")
jpeg(filename = "M_estimator.jpg",
     width = 6, height = 4, units = "in", res = 600)

par(mar = c(4, 4, 3, 1) + 0.1)

plot(squared.loss,
     xlim = c(-10, 10),
     ylim = c(0, 25),
     xlab = "Residual",
     ylab = "Loss value",
     main = "Robust methods",
     bty = "l",
     lwd = 2)
curve(absolute.loss, add = TRUE, lwd = 2, col = "blue")
curve(huber.loss, add = TRUE, lwd = 2, col = "red", lty = 2)
curve(tukey.biweight.loss, add = TRUE, lwd = 2, col = "darkgreen", lty = 2)
legend(x = "bottomright", legend = c("Squared", "Absolute", "Huber", "Biweight"),
       lty = c(1, 1, 2, 2),
       lwd = c(2, 2, 2, 2),
       cex = 0.8,
       col = c("black", "blue", "red", "darkgreen"),
       bg = "white")

dev.off()
