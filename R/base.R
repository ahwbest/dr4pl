# -----------------------------------------------------------------------------
### Dose response relation analysis (drra)
### Hyowon An, UNC Lineberger Comprehensive Cancer Center
### Last updated: 09/02/2016
#

#' Compute predicted responses.
#
#' @param x Dose
#' @param theta Parameters
#
#' @return Predicted response values.
#' @export
MeanResponseCurve <- function(x, theta) {

  f <- theta[1] + (theta[4] - theta[1])/(1 + (x/theta[2])^theta[3])

  return(f)
}

#' Squares of residuals
#'
#' @param r Residuals
#
#' @return Squared residuals
#' @export
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
#' @return result: Huber's loss function values evaluated at residuals r.
#' @export
HuberLoss <- function(r) {
  # This value should be chosen in an adaptive fashion.
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
#' @export
TukeyBiweightLoss <- function(r) {

  const <- 1.345

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
#' @export
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

    n <- length(y)
    f <- theta[1] + (theta[4] - theta[1])/(1 + (x/theta[2])^theta[3])

    return(sum(loss.fcn(y - f))/n)
  }

  return(err.fcn)
}

#' Compute gradient values.
#
#' @param theta Parameters
#' @param dose Dose
#' @param response Response
#
#' @return Gradient values.
#' @export
GradientFunction <- function(theta, dose, response) {
  x <- dose
  y <- response

  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]

  eta <- (x/theta.2)^theta.3
  f <- theta.1 + (theta.4 - theta.1)/(1 + eta)

  deriv.f.theta.1 <- 1 - 1/(1 + eta)
  deriv.f.theta.2 <- (theta.4 - theta.1)*theta.3/theta.2*eta/(1 + eta)^2
  deriv.f.theta.3 <- -(theta.4 - theta.1)/theta.3*log(eta)*eta/(1 + eta)^2
  deriv.f.theta.4 <- 1/(1 + eta)

  deriv.f.theta.2[eta == Inf] <- 0
  deriv.f.theta.3[eta == Inf] <- 0

  return(-2*(y - f)%*%cbind(deriv.f.theta.1, deriv.f.theta.2, deriv.f.theta.3, deriv.f.theta.4))
}

#' Compute the Jacobian matrix
#
#' @param theta Parameters
#' @param x FILL ME OUT!
#'
#' @return Jacobian matrix
#' @export
DerivativeF <- function(theta, x) {

  eta <- (x/theta[2])^theta[3]
  f <- theta[1] + (theta[4] - theta[1])/(1 + eta)

  deriv.f.theta.1 <- 1 - 1/(1 + eta)
  deriv.f.theta.2 <- (theta[4] - theta[1])*theta[3]/theta[2]*eta/(1 + eta)^2
  deriv.f.theta.3 <- -(theta[4] - theta[1])/theta[3]*log(eta)*eta/(1 + eta)^2
  deriv.f.theta.4 <- 1/(1 + eta)

  deriv.f.theta.2[eta == Inf] <- 0
  deriv.f.theta.3[eta == Inf] <- 0

  return(cbind(deriv.f.theta.1, deriv.f.theta.2, deriv.f.theta.3, deriv.f.theta.4))
}

