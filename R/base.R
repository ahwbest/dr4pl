
#library(tensor)

#' Compute predicted responses.
#'
#' @param x Dose
#' @param theta Parameters
#'
#' @return Predicted response values.
#' @export
MeanResponse <- function(x, theta) {

  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]
  
  if(any(is.na(theta))) {
    
    stop("One of the parameter values is NA.")
  }
  
  if(theta.2 < 0) {

    stop("The IC50 parameter estimates become negative during the optimization process.")
  }

  f <- theta.1 + (theta.4 - theta.1)/(1 + (x/theta.2)^theta.3)

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
#' @export
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

#' Returns an error function for given robust fitting method
#
#' @param method.robust NULL, absolute, Huber, or Tukey
#'      - NULL: Sum of squares loss
#'      - absolute: Absolute deviation loss
#'      - Huber: Huber's loss
#'      - Tukey: Tukey's biweight loss
#
#' @return Value of the sum of squared residuals
ErrFcn <- function(method.robust) {

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

  err.fcn <- function(theta, x, y) {

    if(length(theta) != 4) {
      stop("The number of parameters is not 4.")
    }

    if(length(x) != length(y)) {
      stop("The numbers of dose values and response values should be the same.")
    }

    n <- length(y)
    f <- MeanResponse(x, theta)

    if(anyNA(f)) {
      stop("Some of the evaluated function values are NA's.")
    }

    return(sum(loss.fcn(y - f))/n)
  }

  return(err.fcn)
}

#' Compute the derivative values of the mean response function.
#
#' @param theta Parameters
#' @param x Dose
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
  deriv.f.theta.3 <- -(theta.4 - theta.1)/theta.3*log(eta)*eta/(1 + eta)^2
  deriv.f.theta.4 <- 1/(1 + eta)
  
  ### Handle the cases when dose values are zeros
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
  
  return(cbind(deriv.f.theta.1, deriv.f.theta.2, deriv.f.theta.3, deriv.f.theta.4))
}

#' Compute gradient values of the sum-of-squares loss function.
#'
#' @param theta Parameters
#' @param x Dose
#' @param y Response
#'
#' @return Gradient values of the sum-of-squares loss function.
GradientSquaredLoss <- function(theta, x, y) {

  f <- MeanResponse(x, theta)  # Mean response values
  n <- length(x)  # Number of data observations
  
  return(-2*(y - f)%*%DerivativeF(theta, x)/n)
  
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
  if(theta.3 < 0) {
    
    second.deriv.f[, , x == 0] <- 0
    
  }
  
  deriv.f <- DerivativeF(theta, x)
  residuals <- Residual(theta, x, y)
  
  hessian <- (2*t(deriv.f)%*%deriv.f -
              2*tensor(A = second.deriv.f, B = residuals, alongA = 3, alongB = 1))/n

  return(hessian)
  
}

#' Compute residuals.
#' 
#' @param x Doses
#' @param y Responses
#' @param theta Parameters
#' 
#' @return Residuals
Residual <- function(theta, x, y) {

  f <- theta[1] + (theta[4] - theta[1])/(1 + (x/theta[2])^theta[3])
  
  return(y - f)
}