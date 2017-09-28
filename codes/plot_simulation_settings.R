# -------------------------------------------------------------------------------------
# Load libraries
# -------------------------------------------------------------------------------------
library(plyr)
library(ggplot2)

# -------------------------------------------------------------------------------------
# Draw a 4PL model
# -------------------------------------------------------------------------------------
Draw4PLModel <- function(x) {

  z <- 10^x

  return(parms[1] + (parms[4] - parms[1])/(1 + (z/parms[2])^parms[3]))
}

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")
jpeg(filename = "simulation_setting.jpg",
     width = 8, height = 4, units = "in", res = 600)

par(mar = c(6, 4, 2, 2) + 0.1, mfrow = c(1, 2))
parms = c(100, 0.01, -1, 0)
curve(Draw4PLModel, from = -5, to = 1,
      main = "Dose response curve",
      xlab = "", 
      ylab = "Mean response",
      ylim = c(0, 120),
      xaxt = "n",
      yaxt = "n",
      bty = "l",
      lwd = 2,
      col = "red")
legend(x = "bottomleft",
       cex = 0.8,
       lwd = 2,
       lty = c(1, 2, 1, 2, 1, 2),
       col = c("red", "red", "magenta", "magenta", "blue", "blue"),
       legend = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"),
       bty = "n")
axis(side = 1, at = seq(from = -5, to = 1, by = 1),
     labels = 10^seq(from = -5, to = 1, by = 1),
     las = 2)
axis(side = 2, at = seq(from = 0, to = 120, by = 20),
     labels = seq(from = 0, to = 120, by = 20),
     las = 2)
mtext(side = 1, text = "Dose", line = 4)

parms = c(100, 0.01, -2, 0)
curve(Draw4PLModel, from = -5, to = 5,
      add = TRUE, lwd = 2, lty = 2, col = "red")

parms = c(100, 0.1, -1, 0)
curve(Draw4PLModel, from = -5, to = 5,
      add = TRUE, lwd = 2, col = "magenta")
parms = c(100, 0.1, -2, 0)
curve(Draw4PLModel, from = -5, to = 5,
      add = TRUE, lwd = 2, lty = 2, col = "magenta")

parms = c(100, 1, -1, 0)
curve(Draw4PLModel, from = -5, to = 5,
      add = TRUE, lwd = 2, col = "blue")
parms = c(100, 1, -2, 0)
curve(Draw4PLModel, from = -5, to = 5,
      add = TRUE, lwd = 2, lty = 2, col = "blue")

levels.dose <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
parms = c(100, 1, -1, 0)
stdev <- 7.5
doses <- rep(levels.dose, each = 5)
responses <- MeanResponse(doses, parms) + rnorm(length(doses), mean = 0, sd = stdev)

plot(x = log10(doses), y = responses,
     xlim = c(-5, 1),
     ylim = c(0, 120),
     pch = 16,
     main = "Example data (Model 5)",
     xlab = "",
     ylab = "Response",
     bty = "l",
     xaxt = "n",
     yaxt = "n")
curve(Draw4PLModel,
      add = TRUE,
      lwd = 2, col = "blue")
axis(side = 1, at = seq(from = -5, to = 1, by = 1),
     labels = 10^seq(from = -5, to = 1, by = 1),
     las = 2)
axis(side = 2, at = seq(from = 0, to = 120, by = 20),
     labels = seq(from = 0, to = 120, by = 20),
     las = 2)
mtext(side = 1, text = "Dose", line = 4)

dev.off()
