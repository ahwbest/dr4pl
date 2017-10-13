# ----------------------------------------------------------------------------
### Load libraries
#
library(gridExtra)

# ----------------------------------------------------------------------------
### Load data
#
n.error.case <- 4
data.error.list <- vector("list", n.error.case)

data.error.list[[1]] <- drc_error_1
data.error.list[[2]] <- drc_error_2
data.error.list[[3]] <- drc_error_3
data.error.list[[4]] <- drc_error_4

# -----------------------------------------------------------------------------
### Make scatter plots without dose-response curves
#
wd <- getwd()
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")

jpeg(filename = "drc_error_cases.jpg",
     width = 10, height = 10, units = "in", res = 600)

ggplot.list <- vector("list", 4)

for(i in 1:n.error.case) {
  
  data.error <- data.error.list[[i]]
  x <- data.error$Dose
  y <- data.error$Response
  n <- length(y)
  
  dr4pl.error <- dr4pl(Response ~ Dose,
                       data = data.error)
  
  theta <- dr4pl.error$par
  residuals <- Residual(theta, x, y)
  
  # We use the median absolute deviation (mad) as a robust estimator of scale 
  # instead of the estimator suggested in Motulsky and Brown (2006)
  # robust.scale <- quantile(abs(residuals), 0.6827)*n/(n - 4)
  scale.robust <- mad(residuals)  
  
  abs.res.sorted <- sort(abs(residuals), index.return = TRUE)$x
  indices.sorted <- sort(abs(residuals), index.return = TRUE)$ix
  
  Q <- 0.05  # Refer to Motulsky and Brown (2006)
  alphas <- Q*seq(from = n, to = 1, by = -1)/n
  p.values <- 2*pt(q = abs.res.sorted/scale.robust, df = n - 4, lower.tail = FALSE)
  
  indices.FDR <- which(p.values < alphas)
  
  if(length(indices.FDR) == 0) {
    
    indices.outlier <- NULL
    
  } else {
    
    indices.outlier <- indices.sorted[seq(from = min(indices.FDR), to = n, by = 1)]
    
  }
  
  ggplot.list[[i]] <- plot(dr4pl.error,
                           type.curve = "data",
                           text.title = paste("Error Case", i),
                           indices.outlier = indices.outlier)
  
}

grid.arrange(ggplot.list[[1]], ggplot.list[[2]], ggplot.list[[3]], ggplot.list[[4]],
             nrow = 2, ncol = 2)

dev.off()

setwd(wd)

# -----------------------------------------------------------------------------
### Fit a robust regression model and test outliers
#
wd <- getwd()
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")

jpeg(filename = "dr4pl_robust_fit_outlier.jpg",
     width = 10, height = 10, units = "in", res = 600)

ggplot.list <- vector("list", 4)

for(i in 1:n.error.case) {

  data.error <- data.error.list[[i]]
  x <- data.error$Dose
  y <- data.error$Response
  n <- length(y)
  
  dr4pl.robust <- dr4pl(Response ~ Dose,
                        data = data.error,
                        method.robust = "Tukey")
  
  theta <- dr4pl.robust$par
  residuals <- Residual(theta, x, y)
  
  # We use the median absolute deviation (mad) as a robust estimator of scale 
  # instead of the estimator suggested in Motulsky and Brown (2006)
  # robust.scale <- quantile(abs(residuals), 0.6827)*n/(n - 4)
  scale.robust <- mad(residuals)  
  
  abs.res.sorted <- sort(abs(residuals), index.return = TRUE)$x
  indices.sorted <- sort(abs(residuals), index.return = TRUE)$ix
  
  Q <- 0.05  # Refer to Motulsky and Brown (2006)
  alphas <- Q*seq(from = n, to = 1, by = -1)/n
  p.values <- 2*pt(q = abs.res.sorted/scale.robust, df = n - 4, lower.tail = FALSE)
  
  indices.FDR <- which(p.values < alphas)
  
  if(length(indices.FDR) == 0) {
    
    indices.outlier <- NULL
    
  } else {
    
    indices.outlier <- indices.sorted[seq(from = min(indices.FDR), to = n, by = 1)]
    
  }
  
  ggplot.list[[i]] <- plot(dr4pl.robust,
                           text.title = paste("Error Case", i),
                           indices.outlier = indices.outlier)
  
}

grid.arrange(ggplot.list[[1]], ggplot.list[[2]], ggplot.list[[3]], ggplot.list[[4]],
             nrow = 2, ncol = 2)

dev.off()

setwd(wd)
