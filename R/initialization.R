
#' Find initial values for the 4PL model.
#'
#' Args:
#' @param x x values
#' @param y y values
#' @param method.init Initialization method
#' @param method.robust Robust fitting method
#'
#' @return theta.IC50.slope Parameter estimates of the IC50 and slope
#' @export
FindInitialParms <- function(x, y, method.init, method.robust) {

  methods.init <- c("logistic", "Mead")  # Vector of available initialization methods
  
  ### Check whether function arguments are appropriate.
  if(length(x) == 0 || length(y) == 0 || length(x) != length(y)) {

    stop("The same numbers of dose and response values should be supplied.")
  }
  
  if(!is.element(method.init, methods.init)) {
    
    stop("The initialization method name should be one of \'logistic\' and \'Mead\'.")
  }
 
  ### Set the upper and lower asymptotes
  y.range <- diff(range(y))
  
  theta.1.4.zero <- 0.001*y.range  # This value will be added to the maximum and minimum

  y.max <- max(y) + theta.1.4.zero
  y.min <- min(y) - theta.1.4.zero
  
  ### Logistic method
  if(method.init == "logistic") {
    
    theta.1.4.grid.size <- 0.025*y.range  # This value is the size of the grid
    
    grid <- c(-2*theta.1.4.grid.size, -theta.1.4.grid.size, 0, theta.1.4.grid.size,
              2*theta.1.4.grid.size)
    
    theta.1.grid <- y.max + grid
    theta.4.grid <- y.min + grid
    
    theta.mat <- matrix(NA, nrow = length(grid)^2, ncol = 4)
    i.row <- 1

    for(theta.1.init in theta.1.grid) {
      
      for(theta.4.init in theta.4.grid) {
        
        data.hill <- data.frame(x = x, y = y)
        data.hill <- subset(data.hill, subset = (data.hill$y<theta.1.init)&
                                                (data.hill$y>theta.4.init)&
                                                (data.hill$x>0))
        
        # There are not enough data to estimate paramter estimates
        if(nrow(data.hill) <= 4) {
          
          next
        }
        
        data.hill$y.logit <- log((data.hill$y - theta.4.init)/(theta.1.init - data.hill$y))  # Logit transformed responses
        data.hill$x.log <- log(data.hill$x)  # Log transformed doses
        
        lm.hill <- lm(y.logit ~ x.log, data = data.hill)
        
        beta.hat <- lm.hill$coefficients
        theta.3.init <- beta.hat[2]
        theta.2.init <- exp(-beta.hat[1]/theta.3.init)
        
        theta.mat[i.row, ] <- c(theta.1.init, theta.2.init, theta.3.init, theta.4.init)
        i.row <- i.row + 1
      }
    }
   
    theta.mat <- theta.mat[!is.na(theta.mat[, 3])&theta.mat[, 3]<0, ] # This restricts approximation to decline curves only

    if(nrow(theta.mat) == 0) {
      
      stop("The logistic method is not applicable for the input data. Please try
            Mead's method instead.")
    }

    err.fcn <- ErrFcn(method.robust)
    errors <- rep(0, nrow(theta.mat))   # To save the error values
    
    for(i in 1:nrow(theta.mat)) {
      
      theta <- theta.mat[i, ]
      errors[i] <- err.fcn(theta, x, y)
    }
    
    theta.init <- theta.mat[which.min(errors), ]
  
  ### Mead's method
  } else if(method.init == "Mead") {
    
    log.x <- log10(x)
    # y.zero.low <- y
    y.lower.bd <- y.min
    y.zero.low <- y - y.lower.bd

    mu.3.vec <- 10^seq(from = 0.1, to = 2.5, by = 0.1)
    theta.mat <- matrix(0, nrow = length(mu.3.vec), ncol = 4)

    for(i in 1:length(mu.3.vec)) {

      mu.3 <- mu.3.vec[i]

      data.lm <- data.frame(y = y.zero.low,
                            y.gamma.x = y.zero.low*mu.3^log.x,
                            Response = 1)

      lm.init <- stats::lm(Response ~ -1 + y + y.gamma.x, data = data.lm)
      
      lm.init.coef <- lm.init$coefficients
      mu.1 <- lm.init.coef[1]
      mu.2 <- lm.init.coef[2]
      
      if((mu.1 < 0)||(mu.2 < 0)) {  # This restircts approximation to decline curves only
        
        next
      }

      theta.mat[i, 1] <- 1/mu.1
      theta.mat[i, 2] <- (mu.1/mu.2)^(1/log10(mu.3))
      theta.mat[i, 3] <- -log10(mu.3)
    }

    theta.mat[, 1] <- theta.mat[, 1] + y.lower.bd
    theta.mat[, 4] <- theta.mat[, 4] + y.lower.bd

    colnames(theta.mat) <- c("Theta1", "Theta2", "Theta3", "Theta4")
    theta.mat <- theta.mat[theta.mat[, 3] != 0, ] 
    
    if(nrow(theta.mat) == 0) {
      
      stop("Mead's method is not applicable for the inupt data. Please try
           the logistic method instead.")
    }
    
    err.fcn <- ErrFcn(method.robust)
    errors <- rep(0, nrow(theta.mat))   # To save the error values

    for(i in 1:nrow(theta.mat)) {

      theta <- theta.mat[i, ]
      errors[i] <- err.fcn(theta, x, y)
    }

    theta.init <- theta.mat[which.min(errors), ]
  }

  names(theta.init) <- c("theta.1", "theta.2", "theta.3", "theta.4")

  return(theta.init)
}

