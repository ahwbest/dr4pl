# ---------------------------------------------------------------------------
### User-defined functions
#
Draw4PLModel <- function(x) {
  
  z <- 10^x
  
  return(parms[1] + (parms[4] - parms[1])/(1 + (z/parms[2])^parms[3]))
}

# ---------------------------------------------------------------------------
### Load data
#
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

data.error.list <- vector("list", 4)

for(i in 1:4) {
  
  file.name <- paste("drc_error_", i, ".csv", sep = "")
  data.error.list[[i]] <- read.csv(file.name)
}

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")

jpeg(filename = "DR4PL_robust_fit_outlier.jpg",
     width = 8, height = 8, units = "in", res = 600)

par(mfrow = c(2, 2), mar = c(6, 4, 3, 1) + 0.1)

outlier.list <- vector("list", 4)
outlier.list[[1]] <- 102
outlier.list[[2]] <- c(2, 8)
outlier.list[[3]] <- c(90, 101)
outlier.list[[4]] <- c(1, 100)

for(i in 1:4) {
  
  data.error <- data.error.list[[i]]
  
  drra.Mead <- drra(Response ~ Dose,
                    data = data.error,
                    method.init = "Mead",
                    method.robust = "absolute")
  
  col.vec <- rep("black", nrow(data.error))
  col.vec[outlier.list[[i]]] <- "red"
  plot(x = log10(data.error$Dose), y = data.error$Response,
       pch = 16,
       main = paste("Error Case ", i),
       xlab = "",
       ylab = "Response",
       bty = "l",
       xaxt = "n",
       col = col.vec)
  axis(side = 1, at = log10(unique(data.error$Dose)),
       labels = unique(data.error$Dose),
       las = 2)
  mtext(side = 1, text = "Dose", line = 4)
  
  parms <- drra.Mead$parameters
  curve(Draw4PLModel, add = TRUE, col = "blue", lwd = 2)
}

dev.off()