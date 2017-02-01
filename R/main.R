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
                    constrained = constrained,
                    grad,
                    init.parm,
                    method.init,
                    method.optim,
                    method.robust) {

  x <- dose  # Vector of dose values
  y <- response  # Vector of responses

  # Choose the loss function depending on the robust estimation method
  err.fcn <- ErrFcn(method.robust)

  if(!is.null(init.parm)) {  # When initial parameter estimates are given

    theta.init <- init.parm

    # Fit a dose-response model
    if(constrained == FALSE) {

      drr <- optim(par = theta.init,
                   fn = err.fcn,
                   gr = grad,
                   method = method.optim,
                   x = x,
                   y = y,
                   control = list(trace = 1))
    } else {

      constraint.matr <- t(as.matrix(c(0, 0, -1, 0)))

      drr <- constrOptim(theta = theta.init,
                         f = Error,
                         ui = constraint.matr,
                         ci = 0,
                         grad = grad,
                         method = method.optim,
                         dose = x,
                         response = y)
    }

    theta <- drr$par
    error <- drr$value

  } else {  # When initial parameter values are not given

    # Set initial values of parameters
    theta.init <- FindInitialParms(x, y, method.init, method.robust)

    names(theta.init) <- c("Left limit", "IC50", "Slope", "Right limit")

    drr <- optim(par = theta.init,
                 fn = err.fcn,
                 gr = grad,
                 method = method.optim,
                 x = x,
                 y = y,
                 control = list(trace = 1))

    theta <- drr$par
    error <- drr$value
  }

  data.drr <- data.frame(Dose = dose, Response = response)

  list(data = data.drr,
       dose = x,
       response = y,
       parameters = theta,
       error.value = error)
}

drra <- function(x, ...) UseMethod("drra")

#' @describeIn drra Used in the default case, supplying a single dose and response variable
#' @param dose Dose
#' @param response Response
drra.default <- function(dose, response,
                         constrained = FALSE,
                         grad = NULL,
                         init.parm = NULL,
                         method.init = "logistic",
                         method.optim = if(is.null(grad)) "Nelder-Mead" else "BFGS",
                         method.robust = NULL,
                         ...) {

  ### If doses and responses are not numeric, then throw.
  # if(class(dose) != "numeric" || class(response) != "numeric") {
  #   stop("Both doses and responses should be given numeric values.")
  # }

  dose <- as.numeric(dose)
  response <- as.numeric(response)

  drra.obj <- drraEst(dose = dose, response = response,
                      constrained = constrained,
                      grad = grad,
                      init.parm = init.parm,
                      method.init = method.init,
                      method.optim = method.optim,
                      method.robust = method.robust)

  drra.obj$call <- match.call()

  class(drra.obj) <- "drra"
  drra.obj

}

#' @title Fit a 4 parameter logistic (4PL) model to dose-response data.
#'
#' @description A general 4PL model fitting function for analysis of
#'   dose-response relation.
#'
#' @param formula A symbolic description of the model to be fit. Either of the
#'   form 'response ~ dose' or as a data frame with response values in first
#'   column and dose values in second column.
#' @param data A data frame containing variables in the model.
#' @param grad A function returning the gradient values for the optimization
#'   methods "BFGS", "CG" and "L-BFGS-B". If it is NULL, the Nelder-Mead method
#'   will be applied.
#' @param init.parm A vector of initial parameters to be optimized in the model.
#' @param method.init The method of obtaining initial values of the parameters.
#'   If it is NULL, a default "logistic" regression method will be used.
#' @param method.optim The method of optimization of the parameters. This method
#'   name is directly applied to the \code{constrOptim} function provided in the
#'   "base" package of R.
#' @param method.robust The robust estimation method to be used to fit a model.
#'   If it is NULL, then the default least sum of squares estimator is used.
#' @param ... Further arguments to be passed to \code{constrOptim}.
#' @return A 'drra' object for which "confint", "gof", "print" and "summary"
#'   methods are implemented. For details, see the help page of each method.
#'   For example, type \code{?confint.drra} to obtain the confidence intervals
#'   of parameters of the 'drra' object.
#' @details This function fits a 4 parameter logistic (4PL) model to dose-response
#'   data. A formula of the model is
#'   \deqn{\theta[1]+(\theta[4]-\theta[1])/(1+(x/\theta[2])^\theta[3])}
#'
#'   \code{method.init} specifies an initialization method to get initial parameter
#'   estimates based on data. The currently supported initialization methods are
#'   "Logistic", "Mead", "Anke". Details of the methods are given in Section
#'   "Details".
#'
#'   \code{method.optim} specifies an optimization method to be used in
#'   "constrOptim" function. The currently supported optimization techniques
#'   include "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent". For
#'   further details, see the help page of \code{\link[stats]{optim}}.
#'
#'   \code{method.robust} chooses a robust estimation method among 4 methods.
#'   The method of estimation is usually identified by the loss function of the
#'   method. This package implements 4 loss functions: sum of squares loss,
#'   absolute deviation loss, Huber's loss and Tukey's biweight loss. Each of
#'   loss function is explained in detail in the vignette.
#' @examples
#' ryegrass.drra <- drra(rootl ~ conc, data = ryegrass)
#'
#' ryegrass.drra
#' @author Hyowon An, Dirk P. Dittmer and J. S. Marron
#' @seealso \code{\link{confint.drra}}, \code{\link{gof.drra}},
#' \code{\link{print.drra}}, \code{\link{summary.drra}}
drra.formula <- function(formula,
                         constrained = FALSE,
                         data = list(),
                         grad = NULL,
                         init.parm = NULL,
                         method.init = "logistic",
                         method.optim = if(is.null(grad)) "Nelder-Mead" else "BFGS",
                         method.robust = NULL,
                         ...) {

  mf <- model.frame(formula = formula, data = data)
  dose <- model.matrix(attr(mf, "terms"), data = mf)[, 2]
  response <- model.response(mf)

  est <- drra.default(dose = dose, response = response,
                      constrained = constrained,
                      grad = grad,
                      init.parm = init.parm,
                      method.init = method.init,
                      method.optim = method.optim,
                      method.robust = method.robust,
                      ...)

  est$call <- match.call()
  est$formula <- formula
  names(est$parameters) <- c("Left limit", "IC50", "Slope", "Right limit")
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

#' @description A default plotting function for a `drra' object.
#'
#' @param object A `drra' object whose mean response function should be plotted.
#' @examples
#' ryegrass.drra <- drra(rootl ~ conc, data = ryegrass)
#'
#' plot(ryegrass.drra)
plot.drra <- function(object, ...) {

  a <- ggplot(aes(x = object$Data$Dose, y = object$Data$Response), data = object$data)

  a <- a + stat_function(fun = MeanResponse,
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

#' Print the drra object to screen.
#'
#' @param x A drra object.
#' @examples
#' ryegrass.drra <- drra(rootl ~ conc,
#'                       data = ryegrass)
#'
#' print(ryegrass.drra)
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
