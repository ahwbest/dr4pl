#' @name dr4pl
#' @docType package
#' @title Dose response data analysis using the 4 Parameter Logistic model
#' @description Main functions related to fitting
#' @author Hyowon An, Lineberger Comprehensive Cancer Center
#' @details Last updated: 09/19/2016
#' @import stats
#' @import graphics
#' @import ggplot2

#' @description Compute the confidence intervals of parameter estimates of a fitted
#'   model.
#' @param object An object of the dr4pl class.
#' @return A matrix of the confidence intervals in which each row represents a
#'   parameter and each column represents the lower and upper bounds of the
#'   confidence intervals of the corresponding parameters.
#' @details This function computes the confidence intervals of the parameters of the
#'   4PL model based on the second order approximation to the Hessian matrix of the
#'   loss function of the model. Refer to Subsection 5.2.2 of 
#'   Seber, G. A. F. and Wild, C. J. (1989). Nonlinear Regression. Wiley Series in
#'   Probability and Mathematical Statistics: Probability and Mathematical
#'   Statistics. John Wiley & Sons, Inc., New York.
#' @examples
#'   obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_1)
#'
#' confint(obj.dr4pl)
#' @author Hyowon An and Justin T. Landis
#' @export
confint.dr4pl <- function(object, ...) {
  
  x <- object$data$Dose
  y <- object$data$Response
  theta <- object$parameters
  hessian <- object$hessian
  
  n <- length(y)  # Number of observations in data
  f <- MeanResponse(x, theta)

  C.hat.inv <- solve(hessian/2)
  s <- sqrt(sum((y - f)^2)/(n - 4))
  
  q.t <- qt(0.975, df = n - 4)
  std.err <- s*sqrt(diag(C.hat.inv))  # Standard error
  ci <- cbind(theta - q.t*std.err, theta + q.t*std.err)
  
  return(ci)
}

dr4plEst <- function(dose, response,
                    constrained = constrained,
                    grad,
                    init.parm,
                    method.init,
                    method.optim,
                    method.robust,
                    trace = 0) {

  convergence <- TRUE
  
  x <- dose  # Vector of dose values
  y <- response  # Vector of responses
  n <- length(x)  # Number of observations

  # Choose the loss function depending on the robust estimation method
  err.fcn <- ErrFcn(method.robust)

  if(!is.null(init.parm)) {  # When initial parameter estimates are given

    # Use given initial parameter estimates
    theta.init <- init.parm
    names(theta.init) <- c("Upper limit", "IC50", "Slope", "Lower limit")

    constr.matr <- matrix(rbind(c(0, 1, 0, 0), c(0, 0, -1, 0)),
                          nrow = 2,
                          ncol = 4)
    constr.vec <- c(0, 0)

    if(any(constr.matr%*%theta.init<constr.vec)) {
      
      stop("Initial parameter values are not in the interior of the feasible region.")
    }
        
    # Fit a dose-response model. The Hill bounds are currently not returned.
    drr <- constrOptim(theta = theta.init,
                       f = err.fcn,
                       grad = grad,
                       ui = constr.matr,
                       ci = constr.vec,
                       method = method.optim,
                       hessian = TRUE,
                       x = x,
                       y = y)
    
    error <- drr$value
    hessian <- drr$hessian
    theta <- drr$par

  } else {  # When initial parameter values are not given.

    # Set initial values of parameters.
    theta.init <- FindInitialParms(x, y, method.init, method.robust)
    names(theta.init) <- c("Upper limit", "IC50", "Slope", "Lower limit")

    constr.matr <- matrix(rbind(c(0, 1, 0, 0), c(0, 0, -1, 0)),
                          nrow = 2,
                          ncol = 4)
    constr.vec <- c(0, 0)

    if(any(constr.matr%*%theta.init<constr.vec)) {
      
      stop("Initial parameter values are not in the interior of the feasible region.")
      
    }
    
    # Fit a dose-response model. The Hill bounds are currently not returned.
    drr <- constrOptim(theta = theta.init,
                       f = err.fcn,
                       grad = grad,
                       ui = constr.matr,
                       ci = constr.vec,
                       method = method.optim,
                       hessian = TRUE,
                       x = x,
                       y = y)
    
    error <- drr$value
    hessian <- drr$hessian
    theta <- drr$par
    
  }

  ### If boundaries are hit.
  if(constr.matr%*%theta == constr.vec) {
    
    convergence <- FALSE
    
    err.fcn <- ErrFcn("absolute")
    
    drr.robust <- constrOptim(theta = theta.init,
                              f = err.fcn,
                              grad = grad,
                              ui = constr.matr,
                              ci = constr.vec,
                              method = method.optim,
                              hessian = TRUE,
                              x = x,
                              y = y)
    
    theta <- drr.robust$par
    residuals <- Residual(theta, x, y)
    robust.scale <- quantile(abs(residuals), 0.6827)*n/(n - 4)
    abs.res.sorted <- sort(abs(residuals))
    
  }
  
  data.drr <- data.frame(Dose = dose, Response = response)

  list(convergence = convergence,
       data = data.drr,
       dose = x,
       response = y,
       parameters = theta,
       error.value = error,
       hessian = hessian)
}

#' @description Fit the 4 parameter logistic model to the data using the function `dr4plEst'
#' @param x Object to dr4pl.
#' @param ... arguments passed to coef
#' @export
dr4pl <- function(x, ...) UseMethod("dr4pl")

#' @describeIn  dr4pl Used in the default case, supplying a single dose and response variable
#' @param dose Dose
#' @param response Response
#' @export
dr4pl.default <- function(dose, response,
                         constrained = FALSE,
                         grad = GradientFunction,
                         init.parm = NULL,
                         method.init = "logistic",
                         method.optim = if(is.null(grad)) "Nelder-Mead" else "BFGS",
                         method.robust = NULL,
                         trace = 0,
                         ...) {

  ### Check whether function arguments are appropriate.
  if(length(dose) == 0 || length(response) == 0 || length(dose) != length(response)) {
    
    stop("The same numbers of dose and response values should be supplied.")
  }
  
  dose <- as.numeric(dose)
  response <- as.numeric(response)

  dr4pl.obj <- dr4plEst(dose = dose, response = response,
                      constrained = constrained,
                      grad = GradientFunction,
                      init.parm = init.parm,
                      method.init = method.init,
                      method.optim = method.optim,
                      method.robust = method.robust,
                      trace = trace)

  dr4pl.obj$call <- match.call()

  class(dr4pl.obj) <- "dr4pl"
  return(dr4pl.obj)
}

#' @title Fit a 4 parameter logistic (4PL) model to dose-response data.
#'
#' @description A general 4PL model fitting function for analysis of
#'   dose-response relation.
#'
#' @param  formula A symbolic description of the model to be fit. Either of the
#'   form 'response ~ dose' or as a data frame with response values in first
#'   column and dose values in second column.
#' @param constrained Boolean for whether or not this function is constrained.
#' @param data A data frame containing variables in the model.
#' @param grad The gradient function that returns gradient values of the loss
#'   function specified by \code{method.robust}. If this parameter is NULL, then
#'   the gradient function for the usual sum-of-squares loss functions will be
#'   used.
#' @param init.parm A vector of initial parameters to be optimized in the model.
#' @param method.init The method of obtaining initial values of the parameters.
#'   If it is NULL, a default "logistic" regression method will be used.
#' @param method.optim The method of optimization of the parameters. This method
#'   name is directly applied to the \code{constrOptim} function provided in the
#'   "base" package of R.
#' @param method.robust The robust estimation method to be used to fit a model.
#'      - NULL: Sum of squares loss
#'      - absolute: Absolute deviation loss
#'      - Huber: Huber's loss
#'      - Tukey: Tukey's biweight loss
#' @param trace Nonnegative number indicating the level of detailedness of output
#'   generated by \code{constrOptim}. As its value increases, more detailed output is
#'   printed.
#' @param ... Further arguments to be passed to \code{constrOptim}.
#' 
#' @return A 'dr4pl' object for which "confint", "gof", "print" and "summary"
#'   methods are implemented. For details, see the help page of each method.
#'   For example, type \code{?confint.dr4pl} to obtain the confidence intervals
#'   of parameters of the 'dr4pl' object.
#'   
#' @details This function fits a 4 parameter logistic (4PL) model to dose-response
#'   data. A formula of the model is
#'   \deqn{\theta[1]+(\theta[4]-\theta[1])/(1+(x/\theta[2])^\theta[3])}
#'
#'   \code{method.init} specifies an initialization method to get initial parameter
#'   estimates based on data. The currently supported initialization methods are
#'   'logistic' and 'Mead'. Details of the methods are given in Section
#'   'Details'.
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
#' ryegrass.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_1)
#'
#' ryegrass.dr4pl
#' @author Hyowon An, Dirk P. Dittmer and J. S. Marron
#' @seealso \code{\link{confint.dr4pl}}, \code{\link{gof.dr4pl}},
#' \code{\link{print.dr4pl}}, \code{\link{summary.dr4pl}}
#' @export
dr4pl.formula <- function(formula,
                         constrained = FALSE,
                         data = list(),
                         grad = GradientFunction,
                         init.parm = NULL,
                         method.init = "logistic",
                         method.optim = if(is.null(grad)) "Nelder-Mead" else "BFGS",
                         method.robust = NULL,
                         trace = 0,
                         ...) {

  mf <- model.frame(formula = formula, data = data)
  dose <- model.matrix(attr(mf, "terms"), data = mf)[, 2]
  response <- model.response(mf)

  est <- dr4pl.default(dose = dose, response = response,
                       constrained = constrained,
                       grad = grad,
                       init.parm = init.parm,
                       method.init = method.init,
                       method.optim = method.optim,
                       method.robust = method.robust,
                       trace = trace,
                       ...)

  est$call <- match.call()
  est$formula <- formula
  names(est$parameters) <- c("Upper limit", "IC50", "Slope", "Lower limit")
  
  return(est)
}

#' @description Coefficient of a `dr4pl' object
#' @title coef
#' @name coef.dr4pl
#' @param object A 'dr4pl' object
#' @param ... arguments passed to coef
#
#' @return A vector of parameters
#' @export
coef.dr4pl <- function(object, ...) {
  
  object$parameters
}

#' @description A default plotting function for a `dr4pl' object.
#' @title plot
#' @name plot.dr4pl
#' @param object A `dr4pl' object whose mean response function should be plotted.
#' @param ... All arguments that can normally be passed to plot.
#' @examples
#' ryegrass.dr4pl <- DR4PL::dr4pl(Response ~ Dose, data = sample_data_1)
#'
#' plot(ryegrass.dr4pl)
#' @export
plot.dr4pl <- function(object, ...) {

  a <- ggplot2::ggplot(aes(x = object$data$Dose, y = object$data$Response), data = object$data)

  a <- a + ggplot2::stat_function(fun = MeanResponse,
                         args = list(theta = object$parameters),
                         size = 1.2)

  a <- a + ggplot2::geom_point(size = I(5), alpha = I(0.8), color = "blue")

  a <- a + ggplot2::labs(title = "Dose response curve",
                         x = "Dose",
                         y = "Response")
  
  # Set parameters for the grids
  a <- a + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 16))
  a <- a + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  a <- a + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  a <- a + ggplot2::scale_x_log10()
  a <- a + ggplot2::theme_bw()
  
  # Set parameters for the titles and text / margin(top, right, bottom, left)
  a <- a + ggplot2::theme(plot.title = ggplot2::element_text(size = 20, margin = ggplot2::margin(0, 0, 10, 0)))
  a <- a + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 16, margin = ggplot2::margin(15, 0, 0, 0)))
  a <- a + ggplot2::theme(axis.title.y = ggplot2::element_text(size = 16, margin = ggplot2::margin(0, 15, 0, 0)))
  a <- a + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16))
  a <- a + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 16))

  return(a)
}

#' Print the dr4pl object to screen.
#'
#' @param object a dr4pl object to be printed
#' @param ... all normally printable arguments
#' @examples
#' ryegrass.dr4pl <- dr4pl(Response ~ Dose,
#'                       data = sample_data_1)
#'
#' print(ryegrass.dr4pl)
print.dr4pl <- function(object, ...) {
  
  cat("Call:\n")
  print(object$call)

  cat("\nCoefficients:\n")
  print(object$parameters)
}

#' Print the dr4pl object summary to screen.
#' @param object a dr4pl object to be summarized
#' @param ... all normally printable arguments
print.summary.dr4pl <- function(object, ...) {
  
  cat("Call:\n")
  print(object$call)
  cat("\n")

  printCoefmat(object$coefficients, P.value = TRUE, has.Pvalue = TRUE)
}

#' Print the dr4pl object summary.
#' @param object a dr4pl object to be summarized
#' @param ... all normal summary arguments
#' @export
summary.dr4pl <- function(object, ...) {

  TAB <- cbind(Estimate = object$parameters,
               StdErr = object$std.err,
               t.value = object$t.value,
               p.value = object$p.value)

  res <- list(call = object$call,
              coefficients = TAB)

  class(res) <- "summary.dr4pl"
  res
}

#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_1
#' @name sample_data_1
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_2
#' @name sample_data_2
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_3
#' @name sample_data_3
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_4
#' @name sample_data_4
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_5
#' @name sample_data_5
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_6
#' @name sample_data_6
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_7
#' @name sample_data_7
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_8
#' @name sample_data_8
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_9
#' @name sample_data_9
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_10
#' @name sample_data_10
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_11
#' @name sample_data_11
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_12
#' @name sample_data_12
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These may or may not have numerical errors in other dose-response curve-packages, but definitly not using these methods.
#' @title sample_data_13
#' @name sample_data_13
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These all have numerical errors in other dose-response curve-packages, but not using these methods.
#' This data set exemplifies the case of a single extreme outlier of one dose measurement. 
#' @title Single High Outlier
#' @name drc_error_1
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These all have numerical errors in other dose-response curve-packages, but not using these methods.
#' This data set exemplifies the case of multiple outliers as well as a small number of observations per dose measurement.
#' @title Multiple High Outliers at Different measurements 
#' @name drc_error_2
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These all have numerical errors in other dose-response curve-packages, but not using these methods.
#' This data set exemplifies the case of multiple outliers at a single dose measurement as well as the support problem.
#' @title Support Problem and Outliers at a Single Dose Level
#' @name drc_error_3
#' @docType data
#' @keywords sample_data
NULL
#' These are a handful of experimentally derived datasets from the wet-laboratory.
#' These all have numerical errors in other dose-response curve-packages, but not using these methods.
#' This data set exemplifies the support problem.
#' @title Support Problem
#' @name drc_error_4
#' @docType data
#' @keywords sample_data
NULL