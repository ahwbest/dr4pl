### Hyowon An, UNC Lineberger Comprehensive Cancer Center
### Last updated: 09/02/2016

# Compute predicted responses.
#
# Args:
#   x: Dose
#   theta: Parameters
#
# Returns:
#   Predicted response values.
MeanResponseCurve <- function(x, theta) {

  f <- theta[1] + (theta[4] - theta[1])/(1 + (x/theta[2])^theta[3])

  return(f)
}

# Squares of residuals
#
# Args:
#   r: Residuals
#
# Returns:
#   Squared residuals
SquaredLoss <- function(r) {

  return(r^2)
}

# Absolute values of residuals.
#
# Args:
#   r: Residuals
#
# Returns:
#   Absolute valued residuals
AbsoluteLoss <- function(r) {

  return(abs(r))
}

# Divide each element of the first vector by the corresponding element of the
# second vector.
#
# Args:
#   x: Vector of numbers to be divided.
#   y: Vector of divisors
#
# Returns:
#   Results of division. If a divisor is infinity or negative infinity, the result
#   is zero.
Divide <- function(x, y) {

  results <- x/y
  results[(!is.finite(y))] <- 0

  return(results)
}

# Values of Huber's loss function evaluated at residuals r.
#
# Args:
#   r: Residuals
#
# Returns:
#   result: Huber's loss function values evaluated at residuals r.
HuberLoss <- function(r) {

  const <- 1.345  # This value should be chosen in an adaptive fashion.

  ret.val <- r^2  # Vector of results
  outer.term <- 2*const*abs(r) - const^2

  outer.idx <- (abs(r) > const)

  ret.val[outer.idx] <- outer.term[outer.idx]

  return(ret.val)
}

# Values of Tukey's biweight loss function evaluated at residuals r.
#
# Args:
#   r: Residuals
#
# Returns:
#   result: Tukey's biweight loss function values evaluated at residuals r.
TukeyBiweightLoss <- function(r) {

  const <- 1.345

  ret.val <- (r^6)/(const^4) - 3*(r^4)/(const^2) + 3*r^2

  ret.val[abs(r) > const] <- const^2

  return(ret.val)
}

# Returns an error function for given robust fitting method
#
# Args:
#   theta: Parameters
#   x: Dose
#   y: Response
#
# Returns:
#   Value of the sum of squared residuals
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

    #cat("Parameters = ", theta, "\n")
    #cat("f = ", f, "\n")

    return(sum(loss.fcn(y - f))/n)
  }

  return(err.fcn)
}

# Compute gradient values.
#
# Args:
#   theta: Parameters
#   dose: Dose
#   response: Response
#
# Returns:
#   Gradient values.

GradientFunction <- function(theta, x, y) {

  # The slope parameter is zero if and only if the left and right asymptotes are
  # the same.
  if(xor(theta[3] == 0, theta[1] == theta[4])) {
    stop("The slope parameter is zero if and only if the left and right asymptotes
         are the same.")
  }

  eta <- (x/theta[2])^theta[3]
  f <- theta[1] + (theta[4] - theta[1])/(1 + eta)

  deriv.f.theta.1 <- 1 - 1/(1 + eta)
  deriv.f.theta.2 <- (theta[4] - theta[1])*theta[3]/theta[2]*eta/(1 + eta)^2
  deriv.f.theta.3 <- -(theta[4] - theta[1])/theta[3]*log(eta)*eta/(1 + eta)^2
  deriv.f.theta.4 <- 1/(1 + eta)

  # The limit of a derivative as x tends to zero depends on the sign of the slope
  # parameter.
  if(theta[3] > 0) {

    deriv.f.theta.1[x == 0] <- 0
    deriv.f.theta.2[x == 0] <- 0
    deriv.f.theta.3[x == 0] <- 0
    deriv.f.theta.4[x == 0] <- 1

  } else if(theta[3] == 0) {

    deriv.f.theta.1[x == 0] <- 0
    deriv.f.theta.2[x == 0] <- (theta[4] - theta[1])*theta[3]/(4*theta[2])
    deriv.f.theta.3[x == 0] <- 0
    deriv.f.theta.4[x == 0] <- 1/2

  } else if(theta[3] < 0) {

    deriv.f.theta.1[x == 0] <- 1
    deriv.f.theta.2[x == 0] <- 0
    deriv.f.theta.3[x == 0] <- 0
    deriv.f.theta.4[x == 0] <- 0

  }

  return(-2*(y - f)%*%cbind(deriv.f.theta.1, deriv.f.theta.2, deriv.f.theta.3, deriv.f.theta.4))
}

# Compute the Jacobian matrix
#
# Args:
#   theta: Parameters
#
# Returns:
#   Jacobian matrix
DerivativeF <- function(theta, x) {

  # The slope parameter is zero if and only if the left and right asymptotes are
  # the same.
  if(xor(theta[3] == 0, theta[1] == theta[4])) {
    stop("The slope parameter is zero if and only if the left and right asymptotes
          are the same.")
  }

  eta <- (x/theta[2])^theta[3]
  f <- theta[1] + (theta[4] - theta[1])/(1 + eta)

  deriv.f.theta.1 <- 1 - 1/(1 + eta)
  deriv.f.theta.2 <- (theta[4] - theta[1])*theta[3]/theta[2]*eta/(1 + eta)^2
  deriv.f.theta.3 <- -(theta[4] - theta[1])/theta[3]*log(eta)*eta/(1 + eta)^2
  deriv.f.theta.4 <- 1/(1 + eta)

  # The limit of a derivative as x tends to zero depends on the sign of the slope
  # parameter.
  if(theta[3] > 0) {

    deriv.f.theta.1[x == 0] <- 0
    deriv.f.theta.2[x == 0] <- 0
    deriv.f.theta.3[x == 0] <- 0
    deriv.f.theta.4[x == 0] <- 1

  } else if(theta[3] == 0) {

    deriv.f.theta.1[x == 0] <- 0
    deriv.f.theta.2[x == 0] <- (theta[4] - theta[1])*theta[3]/(4*theta[2])
    deriv.f.theta.3[x == 0] <- 0
    deriv.f.theta.4[x == 0] <- 1/2

  } else if(theta[3] < 0) {

    deriv.f.theta.1[x == 0] <- 1
    deriv.f.theta.2[x == 0] <- 0
    deriv.f.theta.3[x == 0] <- 0
    deriv.f.theta.4[x == 0] <- 0

  }

  return(cbind(deriv.f.theta.1, deriv.f.theta.2, deriv.f.theta.3, deriv.f.theta.4))
}

# Compute the Hessian matrix
#
# Args:
#    theta: Parameters
#    x: Dose
#
# Returns:
#    A Hessian matrix
Hessian <- function(theta, x, y) {

  n <- length(x)  # Number of observations
  p <- length(theta)  # Number of parameters

  # Second order derivatives of f
  second.deriv.f <- array(dim = c(p, p, n))

  # eta: terms needed in the Hessian matrix computation
  eta <- (x/theta[2])^theta[3]
  eta[x == 0] <- 0

  deriv.eta.2 <- -theta[3]/theta[2]*eta
  deriv.eta.3 <- eta*log(x/theta[2])
  if(theta[3] < 0) {
    deriv.eta.3[x == 0] <- -Inf
  }else {
    deriv.eta.3[x == 0] <- Inf
  }

  second.deriv.f[1, 1, ] <- 0
  second.deriv.f[1, 2, ] <- deriv.eta.2/(1+eta)^2
  second.deriv.f[1, 3, ] <- deriv.eta.3/(1+eta)^2
  second.deriv.f[1, 4, ] <- 0

  second.deriv.f[2, 1, ] <- -theta[3]/theta[2]*eta/(1 + eta)^2
  second.deriv.f[2, 2, ] <- (theta[4] - theta[1])*theta[3]/theta[2]/(1 + eta)^2*
    (-eta/theta[2] + (1 - eta)/(1 + eta)*deriv.eta.2)
  second.deriv.f[2, 3, ] <- (theta[4] - theta[1])/theta[2]/(1 + eta)^2*
    (eta + theta[3]*(1 - eta)/(1 + eta)*deriv.eta.3)
  second.deriv.f[2, 4, ] <- theta[3]/theta[2]*eta/(1 + eta)^2

  second.deriv.f[3, 1, ] <- log(x/theta[2])*eta/(1 + eta)^2
  second.deriv.f[3, 1, x == 0] <- -Inf
  second.deriv.f[3, 2, ] <- (theta[4] - theta[1])/(1 + eta)^2*
    (eta/theta[2] - (log(x/theta[2]))*(1 - eta)/(1 + eta)*deriv.eta.2)
  if(theta[4] > theta[1]) {
    second.deriv.f[3, 2, x == 0] <- Inf
  }else {
    second.deriv.f[3, 2, x == 0] <- -Inf
  }

  second.deriv.f[3, 3, ] <- -(theta[4] - theta[1])*(log(x/theta[2]))*
    (1 - eta)/(1 + eta)^3*deriv.eta.3
  second.deriv.f[3, 3, x == 0] <- Inf
  second.deriv.f[3, 4, ] <- -log(x/theta[2])*eta/(1 + eta)^2
  second.deriv.f[3, 4, x == 0] <- Inf

  second.deriv.f[4, 1, ] <- 0
  second.deriv.f[4, 2, ] <- -deriv.eta.2/(1+eta)^2
  second.deriv.f[4, 3, ] <- -deriv.eta.3/(1+eta)^2
  second.deriv.f[4, 4, ] <- 0

  deriv.f <- DerivativeF(theta, x)
  residual <- Residual(theta, x, y)

  hessian <- 2*t(deriv.f)%*%deriv.f - 2*tensor(second.deriv.f, residual, 3, 1)

  return(hessian)
}
