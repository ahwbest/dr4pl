
# --------------------------------------------------------------------------------
### Load libraries
#
library(drc)

# --------------------------------------------------------------------------------
### User-defined functions
#
DoseResponseCurve <- function(log10.x) {
  
  x <- 10^log10.x
  
  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]
  
  f <- theta.1 + (theta.4 - theta.1)/(1 + (x/theta.2)^theta.3)
  
  return(f)
}

ErrorFunction <- function(theta, data.drm) {
  x <- data.drm$Conc
  y <- data.drm$Response
  
  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]
  
  f <- theta.1 + (theta.4 - theta.1)/(1 + (x/theta.2)^theta.3)
  
  return(sum((y - f)^2))
}

# --------------------------------------------------------------------------------
### Generate figures
#
n.fig <- 57

for(i in 1:ceiling(n.fig/2)) {
  setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Presentation")
  setEPS()
  postscript(paste("drc_constrOptim_comparison_", i, ".eps", sep = ""),
             width = 8,
             height = 4,
             horizontal = FALSE)

  par(mfrow = c(1, 2))

  ### (2*i - 1)-th data set
  setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")
  file.name <- paste("Dirk_data_", 2*i - 1, ".csv", sep = "")

  data.org <- read.csv(file.name)

  x <- data.org$Conc
  y <- data.org$Response
  
  drm.drc <- drm(Response ~ Dose,
                 data = data.org,
                 fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50")))
  
  slope <- drm.drc$coefficients[1]
  lower.limit <- drm.drc$coefficients[2]
  upper.limit <- drm.drc$coefficients[3]
  ic50 <- drm.drc$coefficients[4]
  
  theta.1.drc <- max(lower.limit, upper.limit)
  theta.2.drc <- ic50
  theta.3.drc <- -abs(slope)
  theta.4.drc <- min(lower.limit, upper.limit)
  
  theta.init.drc <- c(theta.1.drc, theta.2.drc, theta.3.drc, theta.4.drc)
  constraint.matr <- t(as.matrix(c(0, 0, -1, 0)))
  
  drm.cO <- constrOptim(theta = theta.init.drc,
                        f = ErrorFunction, 
                        ui = constraint.matr,
                        ci = 0,
                        data.drm = data.org,
                        method = "Nelder-Mead")
  
  theta.fitted.cO <- drm.cO$par

  # Generate plots
  plot(x = log10(x),
       y = y,
       main = paste(2*i - 1, "-th data set", sep = ""),
       xlab = "",
       ylab = "Response",
       lwd = 2,
       pch = 16,
       type = "p",
       bty = "L",
       xaxt = "n")
  
  theta <- theta.fitted.cO
  curve(DoseResponseCurve, type = "l", add = TRUE, lwd = 2, col = "blue",
        lty = 1)
  
  theta <- c(theta.1.drc, theta.2.drc, theta.3.drc, theta.4.drc)
  curve(DoseResponseCurve, type = "l", add = TRUE, lwd = 2, col = "red",
        lty = 2)
  
  mtext(text = "Dose", side = 1, line = 4)
  x.level <- sort(unique(x))
  at.x <- log10(x.level)
  axis(1, at = at.x, labels = 10^at.x,las=2)
  
  legend(x = "topright", col = c("blue", "red"),
         legend = c("cO", "drc"),
         lwd = 2,
         lty = c(1, 2),
         bty = "n")
  
  ### (2*i)-th data set
  if(2*i > n.fig) {
    dev.off()
    break
  }
  
  setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")
  file.name <- paste("drc_data_", 2*i, ".csv", sep = "")
  
  data.org <- read.csv(file.name)
  
  x <- data.org$Conc
  y <- data.org$Response
  
  drm.drc <- drm(Response ~ Conc,
                 data = data.org,
                 fct = LL.4(names=c("Slope","Lower Limit","Upper Limit","IC50")))
  
  slope <- drm.drc$coefficients[1]
  lower.limit <- drm.drc$coefficients[2]
  upper.limit <- drm.drc$coefficients[3]
  ic50 <- drm.drc$coefficients[4]
  
  theta.1.drc <- max(lower.limit, upper.limit)
  theta.2.drc <- ic50
  theta.3.drc <- -abs(slope)
  theta.4.drc <- min(lower.limit, upper.limit)
  
  theta.init.drc <- c(theta.1.drc, theta.2.drc, theta.3.drc, theta.4.drc)
  constraint.matr <- t(as.matrix(c(0, 0, -1, 0)))
  
  drm.cO <- constrOptim(theta = theta.init.drc,
                        f = ErrorFunction, 
                        ui = constraint.matr,
                        ci = 0,
                        data.drm = data.org,
                        method = "Nelder-Mead")
  
  theta.fitted.cO <- drm.cO$par
  
  # Generate plots
  plot(x = log10(x),
       y = y,
       main = paste(2*i, "-th data set", sep = ""),
       xlab = "",
       ylab = "Response",
       lwd = 2,
       pch = 16,
       type = "p",
       bty = "L",
       xaxt = "n")
  
  theta <- theta.fitted.cO
  curve(DoseResponseCurve, type = "l", add = TRUE, lwd = 2, col = "blue",
        lty = 1)
  
  theta <- c(theta.1.drc, theta.2.drc, theta.3.drc, theta.4.drc)
  curve(DoseResponseCurve, type = "l", add = TRUE, lwd = 2, col = "red",
        lty = 2)
  
  mtext(text = "Dose", side = 1, line = 4)
  x.level <- sort(unique(x))
  at.x <- log10(x.level)
  axis(1, at = at.x, labels = 10^at.x,las=2)
  
  legend(x = "topright", col = c("blue", "red"),
         legend = c("cO", "drc"),
         lwd = 2,
         lty = c(1, 2),
         bty = "n")
  
  dev.off()
}