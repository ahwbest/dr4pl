
# -------------------------------------------------------------------------------------
### Load libraries
#
library(plyr)
library(ggplot2)

# -------------------------------------------------------------------------------------
### Draw a 4PL model
#
Draw4PLModel <- function(log10.x) {

  x <- 10^log10.x

  return(parms[1] + (parms[4] - parms[1])/(1 + (x/parms[2])^parms[3]))
}

wd <- getwd()
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")
jpeg(filename = "4PL_different_parameters.jpg",
     width = 8, height = 4, units = "in", res = 600)

### When the IC50 parameter value changes
par(mar = c(6, 4, 1, 1) + 0.1, mfrow = c(1, 2))
parms = c(100, 0.1, -1, 0)
curve(Draw4PLModel, from = -5, to = 5,
      main = "",
      xlab = "",
      ylab = "Mean response",
      xaxt = "n",
      bty = "l",
      lwd = 2,
      col = "red")
legend(x = "topright",
       lwd = 2,
       lty = c(1, 2, 3),
       col = c("red", "magenta", "blue"),
       legend = c(expression(paste(theta[2], "=0.1")),
                  expression(paste(theta[2], "=1")),
                  expression(paste(theta[2], "=10"))),
       bty = "n")
axis(side = 1, at = seq(from = -5, to = 5, by = 1),
     labels = 10^seq(from = -5, to = 5, by = 1),
     las = 2)
mtext(side = 1, text = "Dose", line = 4)

parms = c(100, 1, -1, 0)
curve(Draw4PLModel, from = -5, to = 5,
      add = TRUE,
      lty = 2,
      lwd = 2,
      col = "magenta")

parms = c(100, 10, -1, 0)
curve(Draw4PLModel, from = -5, to = 5,
      add = TRUE,
      lty = 3,
      lwd = 2,
      col = "blue")

### When the slope parameter value changes
parms = c(100, 1, -0.5, 0)
curve(Draw4PLModel, from = -5, to = 5,
      main = "",
      xlab = "",
      ylab = "Mean response",
      ylim = c(0, 100),
      xaxt = "n",
      bty = "l",
      lwd = 2,
      col = "red")
legend(x = "topright",
       lty = c(1, 2, 3),
       lwd = 2,
       col = c("red", "magenta", "blue"),
       legend = c(expression(paste(theta[3], "=-0.5")),
                  expression(paste(theta[3], "=-1")),
                  expression(paste(theta[3], "=-2"))),
       bty = "n")
axis(side = 1, at = seq(from = -5, to = 5, by = 1),
     labels = 10^seq(from = -5, to = 5, by = 1),
     las = 2)
mtext(side = 1, text = "Dose", line = 4)

parms = c(100, 1, -1, 0)
curve(Draw4PLModel, from = -5, to = 5,
      add = TRUE,
      lty = 2,
      lwd = 2,
      col = "magenta")
parms = c(100, 1, -2, 0)
curve(Draw4PLModel, from = -5, to = 5,
      add = TRUE,
      lty = 3,
      lwd = 2,
      col = "blue")

dev.off()

setwd(wd)