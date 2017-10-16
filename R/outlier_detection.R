
#' Detect outliers by the method of Motulsky and Brown (2006).
#' 
#' @param residuals Vector of residuals from a robust fit.
#' 
#' @return Vector of indices of outliers in the input vector of residuals
#' 
#' @details This function detects outliers from a vector of residuals obtained from
#' a robust fit. The method used here is the same with Motulsky and Brown (2006)
#' except that the median absolute deviation is used instead of the sample quantile
#' based estimator suggested in that paper. Based on the False Discovery Rate (FDR)
#' a set of multiple outliers that have lower FDR's than a threshold are reported.
#' 
#' @author Hyowon An
#' 
#' @references \insertRef{Motulsky2006}{dr4pl}
#' 
#' @export
OutlierDetection <- function(residuals) {
  
  n <- length(residuals)  # Number of observations
  
  # We use the median absolute deviation (mad) as a robust estimator of scale 
  # instead of the estimator suggested in Motulsky and Brown (2006)
  scale.robust <- mad(residuals)  
  
  abs.res.sorted <- sort(abs(residuals), index.return = TRUE)$x
  indices.sorted <- sort(abs(residuals), index.return = TRUE)$ix
  
  Q <- 0.01  # Refer to Motulsky and Brown (2006)
  alphas <- Q*seq(from = n, to = 1, by = -1)/n
  p.values <- 2*pt(q = abs.res.sorted/scale.robust, df = n - 4, lower.tail = FALSE)
  
  indices.FDR <- which(p.values < alphas)
  
  if(length(indices.FDR) == 0) {
    
    indices.outlier <- NULL
  } else {
    
    indices.outlier <- indices.sorted[seq(from = min(indices.FDR), to = n, by = 1)]
  }
  
  return(indices.outlier)
}

