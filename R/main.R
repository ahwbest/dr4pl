# -----------------------------------------------------------------------------
### Dose response relation analysis (drra)
### main.R: Main functions related to fitting
### Hyowon An, Lineberger Comprehensive Cancer Center
### Last updated: 09/19/2016
#

# -----------------------------------------------------------------------------
### Load libraries and source codes
#
library(ggplot2)
library(tensor)
source(".\\R\\base.R")
source(".\\R\\initialization.R")

# -----------------------------------------------------------------------------
### Methods
#
drraEst <- function(dose, response,
                    grad,
                    init.parm,
                    method.init,
                    method.optim,
                    method.robust) {
  # Fit the 4 parameter logistic model to data
  #
  # Args:
  #   dose: Dose,
  #   response: Response,
  #   data: Data
  #   grad: Gradient function
  #   init.parm: Vector of initial parameters
  #   method.init: Initialization method
  #   method.optim: Optimization method
  #   method.robust: Robust estimation method
  #      - NULL: Sum of squares loss
  #      - absolute: Absolute deviation loss
  #      - Huber: Huber's loss
  #      - Tukey: Tukey's biweight loss
  #
  # Returns:
  #   The `drra' object containing parameter estimates
  x <- dose
  y <- response

  data.drra <- data.frame(Dose = x, Response = y)

  ### Fit a 4PL model using `optim' function
  if(!is.null(init.parm)) {  # When initial parameter estimates are given
    theta.init <- init.parm
    constraint.matr <- t(as.matrix(c(0, 0, -1, 0)))

    # Fit a 4PL model using the package constrOptim
    drm.cO <- constrOptim(theta = theta.init,
                          f = Error,
                          ui = constraint.matr,
                          ci = 0,
                          grad = grad,
                          method = method.optim,
                          dose = dose,
                          response = response)

    theta <- drm.cO$par
    error <- drm.cO$value

  } else {

    # Set initial values of parameters
    theta.1.4.init <- FindLeftRightAsymptotes(x, y)

    theta.1.init <- theta.1.4.init[1]
    theta.4.init <- theta.1.4.init[2]

    theta.2.3.init <- FindIC50Slope(x, y, theta.1.4.init, method.init)

    theta.init <- c(theta.1.init, theta.2.3.init, theta.4.init)
    names(theta.init) <- c("Left limit", "IC50", "Slope", "Right limit")

    err.fcn <- ErrFcn(method.robust)

    drm.cO <- optim(par = theta.init,
                    fn = err.fcn,
                    gr = grad,
                    method = method.optim,
                    x = x,
                    y = y,
                    control = list(trace = 2))

    theta <- drm.cO$par
    error <- drm.cO$value
  }

  data.drra <- data.frame(Dose = dose, Response = response)

  list(data = data.drra,
       dose = x,
       response = y,
       parameters = theta,
       error.value = error)
}

# drraEst.cO <- function(dose, response,
#                        grad,
#                        init.parm,
#                        method,
#                        method.robust) {
#   # Fit the 4 parameter logistic model to data
#   #
#   # Args:
#   #   dose: Dose values
#   #   response: Response values
#   #
#   # Returns:
#   #   The `drra' object containing parameter estimates
#   x <- dose
#   y <- response
#
#   data.drra <- data.frame(Dose = dose, Response = response)
#
#   ### Fit a 4PL model using `constrOptim' function
#   if(!is.null(init.parm)) {  # When initial parameter estimates are given
#     theta.init <- init.parm
#     constraint.matr <- t(as.matrix(c(0, 0, -1, 0)))
#
#     # Fit a 4PL model using the package constrOptim
#     drm.cO <- constrOptim(theta = theta.init,
#                           f = Error,
#                           ui = constraint.matr,
#                           ci = 0,
#                           grad = grad,
#                           method = method,
#                           dose = dose,
#                           response = response)
#
#     theta <- drm.cO$par
#     error <- drm.cO$value
#
#   } else {
#
#     # Set initial values of parameters
#     scale.inc <- 0.001
#     y.range <- range(y)
#     len.y.range <- scale.inc * diff(y.range)
#
#     theta.1.init <- max(response) + len.y.range
#     theta.4.init <- min(response) - len.y.range
#
#     y.transf <- log10((response - theta.4.init)/(theta.1.init - response))
#     x.log10 <- log10(dose)
#
#     data.lm <- data.frame(x = x.log10, y = y.transf)
#     data.lm <- data.lm[data.lm$x != -Inf, ]
#
#     lm.init <- lm(y ~ x, data = data.lm)  # Linear model for initial parameter estimates
#     beta.hat <- lm.init$coefficients
#
#     theta.3.init <- beta.hat[2]
#     theta.2.init <- 10^(-beta.hat[1]/theta.3.init)
#
#     theta.init <- c(theta.1.init, theta.2.init, theta.3.init, theta.4.init)
#     constraint.matr <- t(as.matrix(c(0, 0, -1, 0)))
#
#     # Fit the 4PL model using the package constrOptim
#     drm.cO <- constrOptim(theta = theta.init,
#                           f = Error,
#                           ui = constraint.matr,
#                           ci = 0,
#                           grad = grad,
#                           method = method,
#                           dose = dose,
#                           response = response)
#
#     theta <- drm.cO$par
#     error <- drm.cO$value
#   }
#
#   data.drra <- data.frame(Dose = dose, Response = response)
#
#   list(data = data.drra,
#        dose = dose,
#        response = response,
#        parameters = theta,
#        error.value = error)
# }

drra <- function(x, ...) UseMethod("drra")

drra.default <- function(dose, response,
                         grad = NULL,
                         init.parm = NULL,
                         method.init = "logistic",
                         method.optim = if(is.null(grad)) "Nelder-Mead" else "BFGS",
                         method.robust = NULL,
                         ...) {
  # Fit the 4PL model using the function `drraEst'
  #
  # Args:
  #   formula: Formula
  #   data: Data
  #   grad: Gradient function
  #   init.parm: Vector of initial parameters
  #   method.init: Initialization method
  #   method.optim: Optimization method
  #   method.robust: Robust estimation method
  #      - NULL: Sum of squares loss
  #      - absolute: Absolute deviation loss
  #      - Huber: Huber's loss
  #      - Tukey: Tukey's biweight loss
  #
  # Returns:
  #   drra.obj: The object of class `drra'
  dose <- as.numeric(dose)
  response <- as.numeric(response)

  drra.obj <- drraEst(dose = dose, response = response,
                      grad = grad,
                      init.parm = init.parm,
                      method.init = method.init,
                      method.optim = method.optim,
                      method.robust = method.robust)

  drra.obj$call <- match.call()

  class(drra.obj) <- "drra"
  drra.obj
}

drra.formula <- function(formula,
                         data = list(),
                         grad = NULL,
                         init.parm = NULL,
                         method.init = "logistic",
                         method.optim = if(is.null(grad)) "Nelder-Mead" else "BFGS",
                         method.robust = NULL,
                         ...) {
  # Fit the 4PL model using the function `drraEst' and a formula
  #
  # Args:
  #   formula: Formula
  #   data: Data
  #   grad: Gradient function
  #   init.parm: Vector of initial parameters
  #   method.init: Initialization method
  #   method.optim: Optimization method
  #   method.robust: Robust estimation method
  #      - NULL: Sum of squares loss
  #      - absolute: Absolute deviation loss
  #      - Huber: Huber's loss
  #      - Tukey: Tukey's biweight loss
  #
  # Returns:
  #   drra.obj: The object of class `drra'
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)[, 2]
  y <- model.response(mf)

  est <- drra.default(dose = x, response = y,
                      grad = grad,
                      init.parm = init.parm,
                      method.init = method.init,
                      method.optim = method.optim,
                      method.robust = method.robust,
                      ...)

  est$call <- match.call()
  est$formula <- formula
  names(est$parameters) <- c("Upper limit", "IC50", "Slope", "Lower limit")
  est
}

coef.drra <- function(object, ...) {
  # Coefficient of a `drra' object
  #
  # Args:
  #   object: A `drra' object
  #
  # Returns:
  #   A vector of parameters
  object$parameters
}

confint.drra <- function(object, ...) {
  x <- object$data$Dose
  y <- object$data$Response
  theta <- object$parameters

  n <- length(y)  # Number of observations in data
  f <- theta[1] + (theta[4] - theta[1])/(1 + (x/theta[2])^theta[3])
  jacobian <- DerivativeF(theta, x)  # Jacobian matrix

  C.hat.inv <- solve(Hessian(theta, x, y)/2)
  s <- sqrt(sum((y - f)^2)/(n - 4))

  q.t <- qt(0.975, df = n - 4)
  std.err <- s*sqrt(diag(C.hat.inv))  # Standard error
  ci <- cbind(theta - q.t*std.err, theta + q.t*std.err)

  return(ci)
}

plot.drra <- function(object, ...) {
  # Make a scatter plot of a `drra' object
  #
  # Args:
  #   object: A `drra' object
  a <- ggplot(aes(x = Dose, y = Response), data = object$data)

  a <- a + stat_function(fun = MeanResponseCurve,
                         args = list(theta = object$parameters),
                         size = 1.2)

  a <- a + geom_point(size = I(5), alpha = I(0.8), color = "blue")

  a <- a + labs(title = "Dose response curve")

  # Set parameters for the grids
  a <- a + theme(strip.text.x = element_text(size = 16))
  a <- a + theme(panel.grid.minor = element_blank())
  a <- a + theme(panel.grid.major = element_blank())
  a <- a + scale_x_log10()
  a <- a + theme_bw()

  # Set parameters for the titles and text / margin(top, right, bottom, left)
  a <- a + theme(plot.title = element_text(size = 20, margin = margin(0, 0, 10, 0)))
  a <- a + theme(axis.title.x = element_text(size = 16, margin = margin(15, 0, 0, 0)))
  a <- a + theme(axis.title.y = element_text(size = 16, margin = margin(0, 15, 0, 0)))
  a <- a + theme(axis.text.x = element_text(size = 16))
  a <- a + theme(axis.text.y = element_text(size = 16))

  plot(a)
}

print.drra <- function(x, ...) {
  # Print a `drra' object to screen
  #
  # Args:
  #   x: A `drra' object
  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients:\n")
  print(x$parameters)
}

print.summary.drra <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")

  printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE)
}

summary.drra <- function(object, ...) {
  # Summary of a `drra' object
  #
  # Args:
  #   object: A `drra' object
  #
  # Returns:
  #   res: A `summary.drra' object
  TAB <- cbind(Estimate = object$parameters,
               StdErr = object$std.err,
               t.value = object$t.value,
               p.value = object$p.value)

  res <- list(call = object$call,
              coefficients = TAB)

  class(res) <- "summary.drra"
  res
}
