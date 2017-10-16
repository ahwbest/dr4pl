
#' Transform parameters of a 4PL model (theta) into re-parameterized parameters
#' (retheta).
#' 
#' @param theta Parameters of a 4PL model in the dose scale.
#' 
#' @return Reparameterized parameters of a 4PL model among which the EC50 parameter
#' is in the log 10 dose scale.
ParmToLog <- function(theta) {

  # Check whether function arguments are appropriate. 
  if(theta[2]<=0) {
    
    stop("The EC50 or IC50 parameter should always be poistive.")
  }
  
  retheta <- c(theta[1], log10(theta[2]), theta[3], theta[4])
  if(theta[3]<=0) {
    
    names(retheta) <- c("UpperLimit", "Log(IC50)", "Slope", "LowerLimit")
  } else {
    
    names(retheta) <- c("UpperLimit", "Log(EC50)", "Slope", "LowerLimit")
  }
  
  return(retheta)
}

#' Transform reparameterized parameters (retheta) back into original parameters 
#' (theta)
#' 
#' @param retheta Parameters of a 4PL model among which the EC50 or IC50 parameter
#' is in the log 10 dose scale.
#' 
#' @return Parameters of a 4PL model among which the EC50 or IC50 parameter is in 
#' the dose scale.
LogToParm <- function(retheta) {
  
  theta <- c(retheta[1], 10^(retheta[2]), retheta[3], retheta[4])
  
  if(theta[3]<=0) {
    
    names(theta) <- c("UpperLimit", "IC50", "Slope", "LowerLimit")
  } else {
    
    names(theta) <- c("UpperLimit", "EC50", "Slope", "LowerLimit")
  }
  
  return(theta)
}

#' Compute an estimated mean response
#'
#' @param x Dose levels
#' @param theta Parameters of the 4PL model
#'
#' @return Predicted response values.
#' @export
MeanResponse <- function(x, theta) {

  ### Check whether function arguments are appropriate
  if(any(is.na(theta))) {
    
    stop("One of the parameter values is NA.")
  }
  if(theta[2]<=0) {
    
    stop("An IC50 estimate should always be positive.")
  }

  f <- theta[1] + (theta[4] - theta[1])/(1 + (x/theta[2])^theta[3])

  return(f)
}

#' Compute an estimated mean response with the logarithm IC50 parameter
#'
#' @param x Dose levels
#' @param retheta Parameters among which the IC50 parameter is logarithmically
#'   transformed
#'
#' @return Predicted response values.
#' @export
MeanResponseLogIC50 <- function(x, retheta) {
  
  ### Check whether function arguments are appropriate
  if(any(is.na(retheta))) {
    
    stop("One of the parameter values is NA.")
  }
  
  f <- retheta[1] + (retheta[4] - retheta[1])/
                     (1 + 10^(retheta[3]*(log10(x) - retheta[2])))
  
  return(f)
}

#' Squares of residuals
#'
#' @param r Residuals
#
#' @return Squared residuals
SquaredLoss <- function(r) {
  return(r^2)
}

#' Absolute values of residuals
#
#' @param r Residuals
#
#' @return Absolute valued residuals
AbsoluteLoss <- function(r) {

  return(abs(r))
}

#' Values of Huber's loss function evaluated at residuals r.
#
#' @param r Residuals
#' @return Huber's loss function values evaluated at residuals r.
HuberLoss <- function(r) {

  # The value 1.345 was suggested by Huber (1964).
  # See Huber, P. J. (1964). Robust Estimation of a Location Parameter. Annals of Statistics 53(1)
  const <- 1.345

  ret.val <- r^2  # Vector of results
  outer.term <- 2*const*abs(r) - const^2

  outer.idx <- (abs(r) > const)

  ret.val[outer.idx] <- outer.term[outer.idx]

  return(ret.val)
}

#' Values of Tukey's biweight loss function evaluated at residuals r.
#
#' @param r Residuals
#
#' @return result: Tukey's biweight loss function values evaluated at residuals r.
TukeyBiweightLoss <- function(r) {

  # The value 4.685 was suggested by Tukey.
  const <- 4.685

  ret.val <- (r^6)/(const^4) - 3*(r^4)/(const^2) + 3*r^2

  ret.val[abs(r) > const] <- const^2

  return(ret.val)
}

#' Returns an loss function for given robust fitting method
#
#' @param method.robust NULL, absolute, Huber, or Tukey
#'      - NULL: Sum of squares loss
#'      - absolute: Absolute deviation loss
#'      - Huber: Huber's loss
#'      - Tukey: Tukey's biweight loss
#
#' @return Value of the sum of squared residuals
ErrFcn <- function(method.robust) {

  ### Check whether function arguments are appropriate.
  if(!(is.null(method.robust)||method.robust == "absolute"||method.robust == "Huber"||
     method.robust == "Tukey")) {
    
    stop("The robust estimation method should be one of NULL, \"absolute\",
         \"Huber\" or \"Tukey\".")
  }
  
  loss.fcn <- c()

  if(is.null(method.robust)) {
    loss.fcn <- SquaredLoss
  } else if(method.robust == "absolute") {
    loss.fcn <- AbsoluteLoss
  } else if(method.robust == "Huber") {
    loss.fcn <- HuberLoss
  } else if(method.robust == "Tukey") {
    loss.fcn <- TukeyBiweightLoss
  }

  err.fcn <- function(retheta, x, y) {

    ### Check whether function arguments are appropriate.
    if(length(retheta) != 4) {
      
      stop("The number of parameters is not 4.")
    }
    if(length(x) != length(y)) {
      
      stop("The numbers of dose levels and responses should be the same.")
    }

    n <- length(y)
    f <- MeanResponseLogIC50(x, retheta)

    if(anyNA(f)) {
      
      stop("Some of the evaluated function values are NA's.")
    }

    return(sum(loss.fcn(y - f))/n)
  }

  return(err.fcn)
}

#' Compute the derivative values of the mean response function.
#
#' @param theta Parameters of the 4PL model
#' @param x Dose levels
#
#' @return Derivative values of the mean response function.
DerivativeF <- function(theta, x) {
  
  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]
  
  eta <- (x/theta.2)^theta.3
  
  ### Compute derivatives
  deriv.f.theta.1 <- 1 - 1/(1 + eta)
  deriv.f.theta.2 <- (theta.4 - theta.1)*theta.3/theta.2*eta/(1 + eta)^2
  deriv.f.theta.3 <- -(theta.4 - theta.1)*log(x/theta.2)*eta/(1 + eta)^2
  deriv.f.theta.4 <- 1/(1 + eta)
  
  # Handle the cases when dose values are zeros
  if(theta.3 > 0) {
    
    deriv.f.theta.1[x == 0] <- 0
    deriv.f.theta.2[x == 0] <- 0
    deriv.f.theta.3[x == 0] <- 0
    deriv.f.theta.4[x == 0] <- 1
  } else if(theta.3 == 0) {
    
    deriv.f.theta.1[x == 0] <- 0
    deriv.f.theta.2[x == 0] <- (theta.4 - theta.1)*theta.3/(4*theta.2)
    deriv.f.theta.3[x == 0] <- 0
    deriv.f.theta.4[x == 0] <- 1/2
  } else if(theta.3 < 0) {
    
    deriv.f.theta.1[x == 0] <- 1
    deriv.f.theta.2[x == 0] <- 0
    deriv.f.theta.3[x == 0] <- 0
    deriv.f.theta.4[x == 0] <- 0
  }
  
  deriv.f.theta <- cbind(deriv.f.theta.1, deriv.f.theta.2, deriv.f.theta.3, deriv.f.theta.4)
  
  # Check whether return values are appropriate
  if(anyNA(deriv.f.theta)) {
    
    stop("Some of the derivative values are NA's.")
  }
  
  return(deriv.f.theta)
}

#' Compute the derivative values of the mean response function with respect to
#'   reparametrized parameters including the log 10 IC50 parameter.
#'
#' @param retheta Parameters obtained from the original parameters by
#'   \code{retheta[2] <- log10(theta[2])}
#' @param x Dose levels
#'
#' @return Derivative values of the mean response function.
DerivativeFLogIC50 <- function(retheta, x) {

  theta <- retheta
  theta[2] <- 10^retheta[2]
  
  deriv.f.theta <- DerivativeF(theta, x)
  deriv.f.retheta <- deriv.f.theta
  deriv.f.retheta[, 2] <- log(10)*theta[2]*deriv.f.theta[, 2]
  
  # Check whether return values are appropriate
  if(anyNA(deriv.f.retheta)) {
    
    stop("Some of the derivative values are NA's.")
  }
  
  return(deriv.f.retheta)
}

#' Compute gradient values of the sum-of-squares loss function.
#'
#' @param retheta Parameters among which the IC50 parameter is logarithmically
#'   transformed
#' @param x Dose
#' @param y Response
#'
#' @return Gradient values of the sum-of-squares loss function.
GradientSquaredLossLogIC50 <- function(retheta, x, y) {

  f <- MeanResponseLogIC50(x, retheta)  # Mean response values
  n <- length(x)  # Number of data observations
  
  return(-2*(y - f)%*%DerivativeFLogIC50(retheta, x)/n)
}

#' Compute the Hessian matrix of the sum-of-squares loss function
#' 
#' @param theta Parameters
#' @param x Doses
#' @param y Response
#'
#' @return Hessian matrix of the sum-of-squares loss function.
#' @export
Hessian <- function(theta, x, y) {

  n <- length(x)  # Number of observations
  p <- length(theta)  # Number of parameters
  
  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]

  # Second order derivatives of f
  second.deriv.f <- array(data = 0, dim = c(p, p, n))

  eta <- (x/theta.2)^theta.3  # Term needed in the Hessian matrix computation

  deriv.eta.2 <- -theta.3/theta.2*eta
  deriv.eta.3 <- eta*log(x/theta.2)
  
  second.deriv.f[1, 1, ] <- 0
  second.deriv.f[1, 2, ] <- deriv.eta.2/(1+eta)^2
  second.deriv.f[1, 3, ] <- deriv.eta.3/(1+eta)^2
  second.deriv.f[1, 4, ] <- 0

  second.deriv.f[2, 2, ] <- (theta.3*(theta.4 - theta.1))/(theta.2^2*(1 + eta)^3)*
                            (theta.2*(1 - eta)*deriv.eta.2 - eta*(1 + eta))
  second.deriv.f[2, 3, ] <- (theta.4 - theta.1)/(theta.2*(1 + eta)^3)*
                            (eta*(1 + eta) + theta.3*(1 - eta)*deriv.eta.3)
  second.deriv.f[2, 4, ] <- theta.3/theta.2*eta/(1 + eta)^2

  second.deriv.f[3, 3, ] <- (theta.4 - theta.1)/(theta.3^2*(1 + eta)^3)*
                            (eta*(1 + eta)*log(eta) - theta.3*(1 + eta) -
                             eta*(1 - eta)*log(eta)^2)
  second.deriv.f[3, 4, ] <- -log(x/theta[2])*eta/(1 + eta)^2

  second.deriv.f[4, 4, ] <- 0

  second.deriv.f <- (second.deriv.f + aperm(second.deriv.f, c(2, 1, 3)))/2
  
  # Substitue limits for second derivatives when dose levels are zero
  if(theta.3 <= 0) {
    
    second.deriv.f[, , x == 0|eta == Inf] <- 0
  } else if(theta.3 > 0) {
    
    second.deriv.f[, , x == 0|eta == 0] <- 0
    second.deriv.f[1, 3, x == 0|eta == 0] <- 0
    second.deriv.f[3, 3, x == 0|eta == 0] <- -(theta.4 - theta.1)/theta.3
  }
  
  deriv.f <- DerivativeF(theta, x)
  residuals <- Residual(theta, x, y)
  
  hessian <- 2*t(deriv.f)%*%deriv.f -
             2*tensor(A = second.deriv.f, B = residuals, alongA = 3, alongB = 1)

  return(hessian)
}

#' Compute the Hessian matrix of the sum-of-squares loss function with
#' reparameterization.
#' 
#' @param retheta Parameters of a 4PL model among which the EC50 parameter is
#' in the log 10 dose scale
#' @param x Doses
#' @param y Response
#'
#' @return Hessian matrix of the sum-of-squares loss function in terms of
#' reparameterized parameters.
#' @export
HessianLogIC50 <- function(retheta, x, y) {
  
  theta <- LogToParm(retheta)  # Original parameters
  
  hessian <- Hessian(theta, x, y)  # Hessian matrix of original parameters
  hessian.re <- hessian  # Hessian matrix of transformed parameters
  hessian.re[2, ] <- log(10)*theta[2]*hessian.re[2, ]
  hessian.re[, 2] <- log(10)*theta[2]*hessian.re[, 2]

  row.names(hessian.re) <- c("Retheta1", "Retheta2", "Retheta3", "Retheta4")
  colnames(hessian.re) <- c("Retheta1", "Retheta2", "Retheta3", "Retheta4")
  
  return(hessian.re)
}

#' Compute residuals.
#' 
#' @param theta Parameters of a 4PL model
#' @param x Vector of doses
#' @param y Vector of responses
#' 
#' @return Vector of residuals
Residual <- function(theta, x, y) {
  
  return(y - MeanResponse(x, theta))
}

#' Compute residuals with a reparameterized model.
#' 
#' @param retheta Parameters of a 4pl model among which the IC50 parameter is
#' transformed into a log 10 scale.
#' @param x Vector of doses
#' @param y Vector of responses
#' 
#' @return Vector of residuals
ResidualLogIC50 <- function(retheta, x, y) {
  
  
  return(y - MeanResponseLogIC50(x, retheta))
}