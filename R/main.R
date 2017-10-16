
#' @name dr4pl
#' 
#' @docType package
#' 
#' @import graphics 
#' @import stats
#' @import ggplot2
#' @import tensor
#' @import Rdpack
NULL
#' @title Fitting 4 Parameter Logistic (4PL) models to dose-response data.
#' 
#' @description This function fits a 4PL model to dose-response data. Users can
#' obtain fitted parameter estimates as return values. Using auxiliary functions
#' provided by this R package, users can plot a fitted dose-response curve and
#' obtain confidence intervals of true parameters. In addition, the goodness-of-fit
#' test for model adequacy of the 4PL models can be performed when replicates are
#' available for each dose level.
#' 
#' @export
dr4pl <- function(...) UseMethod("dr4pl")

#' @describeIn dr4pl General 4PL model fitting function for analysis of
#'   dose-response relation.
#'
#' @param  formula Symbolic description of the model to be fit. Either of the
#' form 'response ~ dose' or as a data frame with response values in first
#' column and dose values in second column.
#' @param data Data frame containing variables in the model.
#' @param init.parm Vector of initial parameters to be optimized in the model.
#' @param trend Indicator of whether a dose-response curve is a decreasing 
#' \eqn{\theta[3]<0} or increasing curve \eqn{\theta[3]>0}. The default is "auto" 
#' which indicates that the trend of the curve is automatically determined by
#' data. The option "decreasing" will impose a restriction \eqn{\theta[3]<=0} 
#' while the option "increasing" will impose a restriction \eqn{\theta[3]>=0} in an 
#' optimization process.
#' @param method.init Method of obtaining initial values of the parameters.
#' If it is NULL, a default "Mead" method will be used. Assign
#' "logistic" to use the logistic method.
#' @param method.optim Method of optimization of the loss function specified by
#' \code{method.robust}. This function argument is directly passed to the function
#' \code{\link[stats]{constrOptim}} which is provided in the \pkg{base} package of R.
#' @param method.robust Parameter to select loss function for the robust estimation 
#' method to be used to fit a model. The argument NULL indicates the sum of squares
#' loss, "absolute" indicates the absolute deviation loss, "Huber" indicates Huber's
#' loss and "Tukey" indicates Tukey's biweight loss.
#' @param use.Hessian Indicator of whether the Hessian matrix (TRUE) or the
#' gradient vector is used in the Hill bounds.
#' @param level Confidence level to be used in Hill bounds computation.
#' @param failure.message Indicator of whether a message indicating attainment of
#' the Hill bounds and possible resolutions will be printed to the console (TRUE)
#' or hidden (FALSE).
#' @param ... Further arguments to be passed to \code{constrOptim}.
#' 
#' @return A 'dr4pl' object for which "confint", "gof", "print" and "summary"
#'   methods are implemented. For details, see the help page of each method.
#'   For example, type \code{?confint.dr4pl} to obtain the confidence intervals
#'   of parameters of the 'dr4pl' object.
#'   
#' @details This function fits a 4 parameter logistic (4PL) model to dose-response
#'   data. A formula of the model is
#'   \deqn{\theta[1]+(\theta[4]-\theta[1])/(1+(z/\theta[2])^\theta[3])}
#'
#'   \code{method.init} specifies an initialization method to get initial parameter
#'   estimates based on data. The currently supported initialization methods are
#'   "logistic" and 'Mead'. For further details, see the vignette.
#'
#'   \code{method.optim} specifies an optimization method to be used in
#'   "constrOptim" function. The currently supported optimization techniques
#'   include "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent". For
#'   further details, see the help page of \code{\link[stats]{optim}}.
#'
#'   \code{method.robust} chooses a robust estimation method among 4 methods.
#'   The method of estimation is usually identified by the loss function of the
#'   method. This package supports 4 types of loss functions: sum-of-squares loss,
#'   absolute deviation loss, Huber's loss and Tukey's biweight loss. Each of
#'   loss function is explained in detail in the vignette.
#'   
#' @author Hyowon An, \email{ahwbest@gmail.com}
#' @author Justin T. Landis, \email{jtlandis314@gmail.com}
#' @author Aubrey G. Bailey, \email{aubreybailey@gmail.com}
#' @seealso \code{\link{confint.dr4pl}}, \code{\link{gof.dr4pl}},
#' \code{\link{print.dr4pl}}, \code{\link{summary.dr4pl}}
#' @export
dr4pl.formula <- function(formula,
                          data = list(),
                          init.parm = NULL,
                          trend = "auto",
                          method.init = "Mead",
                          method.robust = NULL,
                          method.optim = "Nelder-Mead",
                          use.Hessian = FALSE,
                          level = 0.9999,
                          failure.message = FALSE,
                          ...) {
  
  mf <- model.frame(formula = formula, data = data)
  dose <- model.matrix(attr(mf, "terms"), data = mf)[, 2]
  response <- model.response(mf)
  
  obj <- dr4pl.default(dose = dose,
                       response = response,
                       init.parm = init.parm,
                       trend = trend,
                       method.init = method.init,
                       method.robust = method.robust,
                       method.optim = method.optim,
                       use.Hessian = use.Hessian,
                       level = level,
                       failure.message = failure.message,
                       ...)
  
  obj$call <- match.call()
  obj$formula <- formula

  return(obj)
}

#' @describeIn dr4pl Used in the default case, supplying a single dose and 
#'   response variable
#'   
#' @param dose Vector of dose levels
#' @param response Vector of responses
#'
#' @examples 
#'   a <- dr4pl(dose = sample_data_1$Dose, 
#'                response = sample_data_1$Response, 
#'                method.init = "logistic")
#'   plot(a)
#'
#'   ##Assign method.init = "Mead" to use Mead's method of estimation. 
#'   # Use method.robust to select desired loss function
#'   b <- dr4pl(Response~Dose, 
#'                data = sample_data_4,
#'                method.init = "Mead", 
#'                method.robust = "Tukey" )
#'   plot(b)
#' 
#'   ##compatable with ggplot
#'   library(ggplot2)
#'   c <- dr4pl(Response~Dose, 
#'              data = drc_error_2,
#'              method.init = "Mead", 
#'              method.robust = "absolute" )
#'   d <- plot(c)
#'   d + scale_x_log10(breaks = c(.00135, .0135, .135, 1.35, 13.5))
#' @export
dr4pl.default <- function(dose,
                          response,
                          init.parm = NULL,
                          trend = "auto",
                          method.init = "Mead",
                          method.robust = NULL,
                          method.optim = "Nelder-Mead",
                          use.Hessian = FALSE,
                          level = 0.9999,
                          failure.message = FALSE,
                          ...) {

  types.trend <- c("auto", "decreasing", "increasing")
  types.method.init <- c("logistic", "Mead")
  types.method.optim <- c("Nelder-Mead", "BFGS", "CG", "SANN")
  
  ### Check errors in functions arguments
  if(!is.numeric(dose)||!is.numeric(response)) {
    
    stop("Both doses and responses should be numeric.")
  }
  if(any(dose<0)) {
    
    stop("Dose levels should be nonnegative.")
  }
  if(length(dose) == 0 || length(response) == 0 || length(dose) != length(response)) {
    
    stop("The same numbers of dose and response values should be supplied.")
  }
  if(!is.element(method.init, types.method.init)) {
    
    stop("The initialization method name should be one of \"logistic\" and \"Mead\".")
  }
  if(!is.element(method.optim, types.method.optim)) {
    
    stop("The optimization method name should be one of \"Nelder-Mead\", \"BFGS\",
         \"CG\", \"L-BFGS-B\" and \"SANN\".")
  }
  if(!is.element(trend, types.trend)) {
    
    stop("The type of the \"trend\" parameter should be one of \"auto\", \"decreasing\" and \"increasing\".")
  }

  # Fit a 4PL model
  obj <- dr4plEst(dose = dose,
                  response = response,
                  init.parm = init.parm,
                  trend = trend,
                  method.init = method.init,
                  method.robust = method.robust,
                  method.optim = method.optim,
                  use.Hessian = use.Hessian,
                  level = level)

  obj$call <- match.call()
  class(obj) <- "dr4pl"
  
  # If any robust estimation method is indicated, report outliers to a user.
  if(!is.null(method.robust)) {
    
    theta <- obj$parameters  # Robust parameter estimates
    residuals <- Residual(theta, dose, response)  # Residuals
    
    indices.outlier <- OutlierDetection(residuals)
    
    obj$idx.outlier <- indices.outlier
    obj$robust.plot <- plot(obj, indices.outlier = indices.outlier)
  }
  
  message.diagnosis <- NULL

  ### When convergence failure happens.
  if(obj$convergence == FALSE) {

    ## Decide the method of robust estimation which is more robust than the method
    ## input by a user.
    if(is.null(method.robust)) {
      
      method.robust.new <- "absolute"
    } else if(is.element(method.robust, c("absolute", "Huber"))) {
      
      method.robust.new <- "Tukey"
    } else {
      
      stop("Convergence failure happened but no resolution could be found.")
    }
    
    n <- obj$sample.size  # Number of data points

    # Fit a 4PL model to data
    obj.robust <- dr4plEst(dose = dose,
                           response = response,
                           init.parm = init.parm,
                           trend = trend,
                           method.init = method.init,
                           method.robust = method.robust.new,
                           method.optim = method.optim,
                           use.Hessian = use.Hessian,
                           level = level)

    obj.robust$call <- match.call()
    class(obj.robust) <- "dr4pl"
    
    ## Detect outliers and report them.
    theta <- obj.robust$parameters  # Robust parameter estimates
    residuals <- Residual(theta, dose, response)  # Residuals
    
    indices.outlier <- OutlierDetection(residuals)
    
    obj.robust$idx.outlier <- indices.outlier
    obj.robust$robust.plot <- plot(obj.robust, indices.outlier = indices.outlier)

    ## Print different messages to the console depending on the convergence success
    ## of a robust fit.
    if(obj.robust$convergence) {
      
      message.diagnosis <- 
      paste("The Hill bounds have been hit during optimization, but other robust ",
            "estimation was succesful.\n",
            "Please refer to \"dr4pl.robust\" variable for diagnosis.\n",
            sep = "")
    } else {
      
      message.diagnosis <- 
      paste("The Hill bounds have been hit during optimization with ",
            obj$method.robust, " and ", obj.robust$method.robust, " methods.\n",                
            "Please try other initialization and robust estimation methods.\n", 
            sep = "")
    }
    
    obj$dr4pl.robust <- obj.robust
    obj$message.diagnosis <- message.diagnosis
  }
  
  if(failure.message&&!is.null(message.diagnosis)) {
    
    cat(message.diagnosis)
  }
  
  return(obj)
}

#' @title Private function to fit the 4PL model to dose-response data
#' 
#' @description Private function that actually fits the 4PL model to data. If the
#'   Hill bounds are attained at the end of optimization processes, then an
#'   indicator of convergence failure so that \code{\link{dr4pl.default}} can
#'   look for a remedy for convergence failure.
#' 
#' @name dr4plEst
#' 
#' @param dose Vector of dose levels
#' @param response Vector of responses
#' @param init.parm Vector of initial parameters of the 4PL model supplied by a
#'   user.
#' @param trend Indicator of whether a dose-response curve is a decreasing 
#' \eqn{\theta[3]<0} or increasing curve \eqn{\theta[3]>0}. The default is "auto" 
#' which indicates that the trend of the curve is automatically determined by
#' data. The option "decreasing" will impose a restriction \eqn{\theta[3]<=0} 
#' while the option "increasing" will impose a restriction \eqn{\theta[3]>=0} in an 
#' optimization process.
#' @param method.init Method of obtaining initial values of the parameters.
#' Should be one of "logistic" for the logistic method or "Mead" for the Mead
#' method. The default option is the Mead method.
#' @param method.robust Parameter to select loss function for the robust estimation 
#' method to be used to fit a model. The argument NULL indicates the sum of squares
#' loss, "absolute" indicates the absolute deviation loss, "Huber" indicates Huber's
#' loss and "Tukey" indicates Tukey's biweight loss.
#' @param method.optim Method of optimization of the parameters. This argument
#' is directly delivered to the \code{constrOptim} function provided in the
#' "base" package of R.
#' @param use.Hessian Indicator of whether the Hessian matrix (TRUE) or the
#' gradient vector is used in the Hill bounds.
#' @param level Confidence level to be used in Hill bounds computation.
#' 
#' @return List of final parameter estimates, name of robust estimation, loss value
#' and so on.
dr4plEst <- function(dose, response,
                     init.parm,
                     trend,
                     method.init,
                     method.optim,
                     method.robust,
                     use.Hessian,
                     level) {
  
  convergence <- TRUE
  x <- dose  # Vector of dose values
  y <- response  # Vector of responses
  n <- length(x)  # Number of observations
  
  # Choose the loss function depending on the robust estimation method
  loss.fcn <- ErrFcn(method.robust)
  # Currently only the gradient function for the squared loss is implemented
  grad <- GradientSquaredLossLogIC50
  
  tuning.barrier <- 1e-04  # Tuning parameter for the log Barrier method
  
  ### When initial parameter estimates are given
  if(!is.null(init.parm)) {

    # Check whether initial parameter estimates satisfy constraints
    if(init.parm[2] <= 0) {
      
      stop("The IC50 parameter should be positive.")
    }
    
    # Use given initial parameter estimates
    retheta.init <- init.parm
    retheta.init[2] <- log10(init.parm[2])
    
    names(retheta.init) <- c("Upper limit", "Log10(IC50)", "Slope", "Lower limit")

    constr.mat <- matrix(c(1, 0, 0, -1), nrow = 1, ncol = 4)
    constr.vec <- 0
    
    # Impose a constraint on the slope parameter based on the function argument
    # "trend".
    if(trend == "decreasing") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, -1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    } else if(trend == "increasing") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, 1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    }
    
    # Fit a 4PL model to data
    optim.dr4pl <- constrOptim(theta = retheta.init,
                               f = loss.fcn,
                               grad = grad,
                               ui = constr.mat,
                               ci = constr.vec,
                               method = method.optim,
                               hessian = TRUE,
                               x = x,
                               y = y)

    loss <- optim.dr4pl$value
    hessian <- optim.dr4pl$hessian
    retheta <- optim.dr4pl$par
    
    theta <- retheta
    theta[2] <- 10^retheta[2]
    
  ### When initial parameter values are not given.
  } else {
    
    ## Obtain initial parameter estimates.
    theta.init <- FindInitialParms(x, y, trend, method.init, method.robust)
    retheta.init <- ParmToLog(theta.init)
    
    Hill.bounds <- FindHillBounds(x, y, theta.init, use.Hessian, level)
    
    constr.mat <- matrix(rbind(c(1, 0, 0, -1),
                               c(0, 1, 0, 0),
                               c(0, -1, 0, 0),
                               c(0, 0, 1, 0),
                               c(0, 0, -1, 0)),
                         nrow = 5,
                         ncol = 4)
    constr.vec <- c(0, Hill.bounds$LogTheta2[1], -Hill.bounds$LogTheta2[2],
                    Hill.bounds$Theta3[1], -Hill.bounds$Theta3[2])

    # Impose a constraint on the slope parameter based on the function argument
    # "trend".
    if(trend == "decreasing") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, -1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    } else if(trend == "increasing") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, 1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    }

    if(any(constr.mat%*%retheta.init<constr.vec)) {
      
      stop("Initial parameter values are not in the interior of the feasible region.")
    }

    # Fit the 4PL model
    optim.dr4pl <- constrOptim(theta = retheta.init,
                               f = loss.fcn,
                               grad = grad,
                               ui = constr.mat,
                               ci = constr.vec,
                               method = method.optim,
                               hessian = TRUE,
                               mu = tuning.barrier,
                               x = x,
                               y = y)
    
    loss <- optim.dr4pl$value
    hessian <- optim.dr4pl$hessian
    retheta <- optim.dr4pl$par
    
    theta <- LogToParm(retheta)
  }
  
  ### If the Hill bounds are hit.
  if(any(abs(constr.mat%*%retheta - constr.vec)<tuning.barrier)) {
    
    convergence <- FALSE
  } 
  
  # Data frame consisting of doses and responses
  data.dr4pl <- data.frame(Dose = dose, Response = response)
  
  name.robust <- method.robust
  if(is.null(method.robust)) {
    
    name.robust <- "squared"
  }
  
  list(convergence = convergence,
       sample.size = n,
       data = data.dr4pl,
       parameters = theta,
       loss.value = loss,
       hessian = hessian,
       method.robust = name.robust)
}