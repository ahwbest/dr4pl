
#' @name dr4pl
#' 
#' @docType package
#' 
#' @import graphics
#' @import ggplot2
#' @import stats
#' @import tensor

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
#' @param decline Indicator of whether the curve is a decline \eqn{\theta[3]<0} 
#' or growth curve \eqn{\theta[3]>0}. The default is "auto" which indicates 
#' that no restriction is imposed on the slope parameter \eqn{\theta[3]}. The
#' option "decline" will impose a restriction \eqn{\theta[3]<=0} while the
#' option "growth" will impose a restriction \eqn{\theta[3]>=0} in an optimization
#' process.
#' @param method.init Method of obtaining initial values of the parameters.
#' If it is NULL, a default "logistic" regression method will be used. Assign
#' "Mead" to use Mead's method.
#' @param method.optim Method of optimization of the loss function specified by
#' \code{method.robust}. This function argument is directly passed to the function
#' \code{\link[stats]{constrOptim}} which is provided in the \pkg{base} package of R.
#' @param method.robust Parameter to select loss function for the robust estimation 
#' method to be used to fit a model. 
#' - NULL: Sum of squares loss 
#' - absolute: Absolute deviation loss 
#' - Huber: Huber's loss 
#' - Tukey: Tukey's biweight loss
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
#'   'logistic' and 'Mead'. For further details, see the vignette.
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
                          decline = "auto",
                          method.init = "logistic",
                          method.optim = "Nelder-Mead",
                          method.robust = NULL,
                          ...) {
  
  mf <- model.frame(formula = formula, data = data)
  dose <- model.matrix(attr(mf, "terms"), data = mf)[, 2]
  response <- model.response(mf)
  
  est <- dr4pl.default(dose = dose,
                       response = response,
                       init.parm = init.parm,
                       decline = decline,
                       method.init = method.init,
                       method.optim = method.optim,
                       method.robust = method.robust,
                       ...)
  
  est$call <- match.call()
  est$formula <- formula
  names(est$parameters) <- c("Upper limit", "IC50", "Slope", "Lower limit")
  
  return(est)
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
                          decline = "auto",
                          method.init = "logistic",
                          method.optim = "Nelder-Mead",
                          method.robust = NULL,
                          ...) {

 
   types.decline <- c("auto", "decline", "growth")
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
  if(!is.element(decline, types.decline)) {
    
    stop("The type of the \"decline\" parameter should be one of \"auto\", \"decline\" and \"growth\".")
  }

  # Fit a 4PL model
  obj.dr4pl <- dr4plEst(dose = dose,
                        response = response,
                        init.parm = init.parm,
                        decline = decline,
                        method.init = method.init,
                        method.optim = method.optim,
                        method.robust = method.robust)

  ### When convergence failure happens.
  if(obj.dr4pl$convergence == FALSE) {  

    
    ### Decide the method of robust estimation which is more robust than the method
    ### input by a user.
    if(is.null(method.robust)) {
      
      method.robust.new <- "absolute"
    } else if(is.element(method.robust, c("absolute", "Huber"))) {
      
      method.robust.new <- "Tukey"
    } else {
      
      stop("Convergence failure happened but no resolution could be found.")
    }
    
    n <- obj.dr4pl$sample.size  # Number of data points
    
    # Fit a 4PL model to data
    obj.dr4pl <- dr4plEst(dose = dose, 
                          response = response,
                          init.parm = init.parm,
                          decline = decline,
                          method.init = method.init,
                          method.optim = method.optim,
                          method.robust = method.robust.new)
    
    theta <- obj.dr4pl$parameters
    residuals <- Residual(theta, dose, response)
    
    # We use the median absolute deviation (mad) as a robust estimator of scale 
    # instead of the estimator suggested in Motulsky and Brown (2006)
    # scale.robust <- quantile(abs(residuals), 0.6827)*n/(n - 4)
    scale.robust <- mad(residuals)  
    
    abs.res.sorted <- sort(abs(residuals), index.return = TRUE)$x
    indices.sorted <- sort(abs(residuals), index.return = TRUE)$ix
    
    Q <- 0.05  # Refer to Motulsky and Brown (2006)
    alphas <- Q*seq(from = n, to = 1, by = -1)/n
    p.values <- 2*pt(q = abs.res.sorted/scale.robust, df = n - 4, lower.tail = FALSE)
    
    indices.FDR <- which(p.values < alphas)
    
    if(length(indices.FDR) == 0) {
      
      indices.outlier <- NULL
    } else {
      
      indices.outlier <- indices.sorted[seq(from = min(indices.FDR), to = n, by = 1)]
    }
    
    plot(obj.dr4pl, indices.outlier = indices.outlier)
  }
  
  obj.dr4pl$call <- match.call()

  class(obj.dr4pl) <- "dr4pl"
  return(obj.dr4pl)
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
#' @param decline Indicator of whether the curve is a decline \eqn{\theta[3]<0} 
#'   or growth curve \eqn{\theta[3]>0}. The default is "auto" which indicates 
#'   that no restriction is imposed on the slope parameter \eqn{\theta[3]}. The
#'   option "decline" will impose a restriction \eqn{\theta[3]<=0} while the
#'   option "growth" will impose a restriction \eqn{\theta[3]>=0} in an optimization
#'   process.
#' @param method.init Method of obtaining initial values of the parameters.
#'   Should be one of "logistic" for the logistic method or "Mead" for the Mead
#'   method. The default option is the logistic method.
#' @param method.optim Method of optimization of the parameters. This argument
#'   is directly delivered to the \code{constrOptim} function provided in the
#'   "base" package of R.
#' @param method.robust Method of robust estimation. Should be one of the followings.
#'      - NULL: Squares loss 
#'      - absolute: Absolute deviation loss 
#'      - Huber: Huber's loss 
#'      - Tukey: Tukey's biweight loss
dr4plEst <- function(dose, response,
                     init.parm,
                     decline,
                     method.init,
                     method.optim,
                     method.robust) {
  
  convergence <- TRUE
  x <- dose  # Vector of dose values
  y <- response  # Vector of responses
  n <- length(x)  # Number of observations
  
  # Choose the loss function depending on the robust estimation method
  err.fcn <- ErrFcn(method.robust)
  # Currently only the gradient function for the squared loss is implemented
  grad <- GradientSquaredLossLogIC50
  
  ### When initial parameter estimates are given
  if(!is.null(init.parm)) {

    # Check whether initial parameter estimates satisfy constraints
    if(init.parm[2] <= 0) {
      
      stop("The IC50 parameter should be positive.")
    }
    
    # Use given initial parameter estimates
    theta.re.init <- init.parm
    theta.re.init[2] <- log10(init.parm[2])
    
    names(theta.re.init) <- c("Upper limit", "Log10(IC50)", "Slope", "Lower limit")
    
    # Impose a constraint on the slope parameter based on the function argument
    # `decline`.
    if(decline == "decline") {
      
      constr.mat <- matrix(c(0, 0, -1, 0), nrow = 1, ncol = 4)
      constr.vec <- 0
    } else if(decline == "growth") {
      
      constr.mat <- matrix(c(0, 0, 1, 0), nrow = 1, ncol = 4)
      constr.vec <- 0
    }
    
    # Fit a 4PL model to data
    if(decline == "auto") {
      
      optim.dr4pl <- optim(par = theta.re.init,
                           fn = err.fcn,
                           gr = GradientSquaredLossLogIC50,
                           method = method.optim,
                           hessian = TRUE,
                           x = x,
                           y = y)
    } else {
      
      optim.dr4pl <- constrOptim(theta = theta.re.init,
                                 f = err.fcn,
                                 grad = grad,
                                 ui = constr.mat,
                                 ci = constr.vec,
                                 method = method.optim,
                                 hessian = TRUE,
                                 x = x,
                                 y = y)
    }
    
    loss <- optim.dr4pl$value
    hessian <- optim.dr4pl$hessian
    theta.re <- optim.dr4pl$par
    
    theta <- theta.re
    theta[2] <- 10^theta.re[2]
    
  ### When initial parameter values are not given.
  } else {
    
    ### Obtain initial values of parameters.
    theta.init <- FindInitialParms(x, y, decline, method.init, method.robust)
    names(theta.init) <- c("Upper limit", "IC50", "Slope", "Lower limit")
    
    theta.re.init <- theta.init
    theta.re.init[2] <- log10(theta.init[2])
    names(theta.re.init) <- c("Upper limit", "Log(IC50)", "Slope", "Lower limit")
    
    ### Compute confidence intervals of the true parameters
    deriv.f <- DerivativeF(theta.init, x)
    residuals <- Residual(theta.init, x, y)
    
    C.hat.inv <- try(solve(t(deriv.f)%*%deriv.f), silent = TRUE)  # Inverse matrix

    if(inherits(C.hat.inv, "try-error")) {
      
      C.hat.Chol <- try(chol(t(deriv.f)%*%deriv.f, silent = TRUE))  # Cholesky decomposition
      
      if(inherits(C.hat.Chol, "try-error")) {
        
        C.hat.Chol <- try(chol(0.99*t(deriv.f)%*%deriv.f + 0.01*diag(dim(deriv.f)[2])))
        
        if(inherits(C.hat.Chol, "try-error")) {
         
           C.hat.Chol <- NULL
        }
      }
      
      if(!is.null(C.hat.Chol)) {
        
        C.hat.inv <- chol2inv(C.hat.Chol)
      } else {
        
        C.hat.inv <- NULL# Proceed with the method of Wang et al. (2010)
      }
    }
    
    s <- sqrt(sum(residuals^2)/(n - 4))
    
    q.t <- qt(0.9999, df = n - 4)
    std.err <- s*sqrt(diag(C.hat.inv))  # Standard error
    ci <- cbind(theta.init - q.t*std.err, theta.init + q.t*std.err)  # Confidence intervals

    ### Perform constrained optimization
    bounds.theta.2 <- ci[2, ]
    bounds.theta.3 <- ci[3, ]
    
    # Sometimes the lower bound of the IC50 parameter is negative.
    if(bounds.theta.2[1]<0) {
      
      constr.mat <- matrix(rbind(c(0, -1, 0, 0),
                                 c(0, 0, 1, 0),
                                 c(0, 0, -1, 0)),
                           nrow = 3,
                           ncol = 4)
      constr.vec <- c(-log10(bounds.theta.2[2]), 
                      bounds.theta.3[1], -bounds.theta.3[2])

    } else {
      
      constr.mat <- matrix(rbind(c(0, 1, 0, 0),
                                 c(0, -1, 0, 0),
                                 c(0, 0, 1, 0),
                                 c(0, 0, -1, 0)),
                           nrow = 4,
                           ncol = 4)
      constr.vec <- c(log10(bounds.theta.2[1]), -log10(bounds.theta.2[2]),
                      bounds.theta.3[1], -bounds.theta.3[2])
    }
    
    # Impose a constraint on the slope parameter based on the function argument
    # `decline`.
    if(decline == "decline") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, -1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    } else if(decline == "growth") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, 1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    }

    if(any(constr.mat%*%theta.re.init<constr.vec)) {
      
      stop("Initial parameter values are not in the interior of the feasible region.")
    }

    # Fit the 4PL model
    optim.dr4pl <- constrOptim(theta = theta.re.init,
                               f = err.fcn,
                               grad = grad,
                               ui = constr.mat,
                               ci = constr.vec,
                               method = method.optim,
                               hessian = TRUE,
                               x = x,
                               y = y)
    
    loss <- optim.dr4pl$value
    hessian <- optim.dr4pl$hessian
    theta.re <- optim.dr4pl$par
    
    theta <- theta.re
    theta[2] <- 10^theta.re[2]
  }
  
  ### If boundaries are hit
  if(any(constr.mat%*%theta.re == constr.vec)) {
    
    convergence <- FALSE
  } 
  
  data.dr4pl <- data.frame(Dose = dose, Response = response)
  
  list(convergence = convergence,
       data = data.dr4pl,
       dose = x,
       response = y,
       sample.size = n,
       parameters = theta,
       loss.value = loss,
       hessian = hessian)
}