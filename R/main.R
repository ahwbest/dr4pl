#' @name drra
#' @docType package
#' @title  Dose response relation analysis (drra)
#' @description Main functions related to fitting
#' @author Hyowon An, Lineberger Comprehensive Cancer Center
#' @details Last updated: 09/19/2016
#' @import stats
#' @import graphics
#' @import ggplot2



#' @param grad Gradient function
#' @param init.parm Vector of initial parameters
#' @param method.init Initialization method
#' @param method.optim Optimization method
#' @param method.robust Robust estimation method
#'      - NULL: Sum of squares loss
#'      - absolute: Absolute deviation loss
#'      - Huber: Huber's loss
#'      - Tukey: Tukey's biweight loss
drraEst <- function(dose, response,
                    grad,
                    init.parm,
                    method.init,
                    method.optim,
                    method.robust) {

  x <- dose
  y <- response

  data.drra <- data.frame(Dose = x, Response = y)

  ### Fit a 4PL model using `optim' function
  if(!is.null(init.parm)) {  # When initial parameter estimates are given
    theta.init <- init.parm
    constraint.matr <- t(as.matrix(c(0, 0, -1, 0)))

    err.fcn <- ErrFcn(method.robust)
    # Fit a 4PL model using the package constrOptim
    drm.cO <- constrOptim(theta = theta.init,
                                 f = err.fcn,
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


#' @description Fit the 4 parameter logistic model to the data using the function `drraEst'
#' @param x Object to drra.
#' @param ... arguments passed to coef
#' @export
drra <- function(x, ...) UseMethod("drra")


#' @describeIn drra Used in the default case, supplying a single dose and response variable
#' @param dose Dose
#' @param response Response
drra.default <- function(dose, response,
                         grad = NULL,
                         init.parm = NULL,
                         method.init = "logistic",
                         method.optim = if(is.null(grad)) "Nelder-Mead" else "BFGS",
                         method.robust = NULL,
                         ...) {

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


#' @describeIn drra Used as a formula call for a list of multiple data.
#' @param data list of dose and response vectors
#' @param formula Formula
#' @export
drra.formula <- function(formula,
                         data = list(),
                         grad = NULL,
                         init.parm = NULL,
                         method.init = "logistic",
                         method.optim = if(is.null(grad)) "Nelder-Mead" else "BFGS",
                         method.robust = NULL,
                         ...) {

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

#' @description Coefficient of a `drra' object
#' @title coef
#' @name coef.drra
#' @param object A 'drra' object
#' @param ... arguments passed to coef
#
#' @return A vector of parameters
#' @export
coef.drra <- function(object, ...) {
  object$parameters
}

plot.drra <- function(object, ...) {

  a <- ggplot(aes(x = object$Data$Dose, y = object$Data$Response), data = object$data)

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

  TAB <- cbind(Estimate = object$parameters,
               StdErr = object$std.err,
               t.value = object$t.value,
               p.value = object$p.value)

  res <- list(call = object$call,
              coefficients = TAB)

  class(res) <- "summary.drra"
  res
}
