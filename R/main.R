#' @name dr4pl
#' @docType package
#' @title Dose-response data analysis using the 4 Parameter Logistic model
#' @description Main functions related to fitting
#' @author Hyowon An, Lineberger Comprehensive Cancer Center
#' @details Last updated: 09/19/2016
#' @import stats
#' @import graphics
#' @import ggplot2
#' @import tensor

# Begin Program

#' @title Fit a 4 parameter logistic (4PL) model to dose-response data.
#' 
#' @export
dr4pl <- function(...) UseMethod("dr4pl")

#' @describeIn dr4pl General 4PL model fitting function for analysis of
#'   dose-response relation.
#'
#' @param  formula Symbolic description of the model to be fit. Either of the
#'   form 'response ~ dose' or as a data frame with response values in first
#'   column and dose values in second column.
#' @param data Data frame containing variables in the model.
#' @param init.parm Vector of initial parameters to be optimized in the model.
#' @param decline Indicator of whether the curve is a decline \eqn{\theta[3]<0} 
#'   or growth curve \eqn{\theta[3]>0}. The default is "auto" which indicates 
#'   that no restriction is imposed on the slope parameter \eqn{\theta[3]}. The
#'   option "decline" will impose a restriction \eqn{\theta[3]<=0} while the
#'   option "growth" will impose a restriction \eqn{\theta[3]>=0} in an optimization
#'   process.
#' @param method.init The method of obtaining initial values of the parameters.
#'   If it is NULL, a default "logistic" regression method will be used. Assign
#'   "Mead" to use Mead's method.
#' @param method.optim The method of optimization of the parameters. This method
#'   name is directly applied to the \code{constrOptim} function provided in the
#'   "base" package of R.
#' @param method.robust Parameter to select error function for the robust estimation method to be used to fit a model. 
#'      - NULL: Sum of squares loss 
#'      - absolute: Absolute deviation loss 
#'      - Huber: Huber's loss 
#'      - Tukey: Tukey's biweight loss
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
#'   
#' @author Hyowon An, Justin T. Landis and Aubrey G. Bailey
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
#'   # Use method.robust to select desired error function
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

  methods.init <- c("logistic", "Mead")
  
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
  if(!is.element(method.init, methods.init)) {
    
    stop("The initialization method name should be one of \'logistic\' and \'Mead\'.")
  }

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
    
    # Fit a dose-response model.
    optim.dr4pl <- optim(par = theta.re.init,
                         fn = err.fcn,
                         gr = grad,
                         method = method.optim,
                         hessian = TRUE,
                         x = x,
                         y = y)

    error <- optim.dr4pl$value
    hessian <- optim.dr4pl$hessian
    theta.re <- optim.dr4pl$par
    
    theta <- theta.re
    theta[2] <- 10^theta.re[2]
    
  ### When initial parameter values are not given.
  } else {
    
    ### Obtain initial values of parameters.
    theta.init <- FindInitialParms(x, y, method.init, method.robust)
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
      constr.vec <- c(-log10(bounds.theta.2[2]), bounds.theta.3[1], -bounds.theta.3[2])

    } else {
      
      constr.mat <- matrix(rbind(c(0, 1, 0, 0),
                                 c(0, -1, 0, 0),
                                 c(0, 0, 1, 0),
                                 c(0, 0, -1, 0)),
                           nrow = 4,
                           ncol = 4)
      constr.vec <- c(log10(bounds.theta.2), bounds.theta.3)
    }

    if(any(constr.mat%*%theta.re.init<constr.vec)) {
      
      stop("Initial parameter values are not in the interior of the feasible region.")
    }

    # Fit the 4PL model
    cO.dr4pl <- constrOptim(theta = theta.re.init,
                            f = err.fcn,
                            grad = grad,
                            ui = constr.mat,
                            ci = constr.vec,
                            method = method.optim,
                            hessian = TRUE,
                            x = x,
                            y = y)
    
    error <- cO.dr4pl$value
    hessian <- cO.dr4pl$hessian
    theta.re <- cO.dr4pl$par
    
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
       error.value = error,
       hessian = hessian)
}

#' @description Default plotting function for a `dr4pl' object. Plot displays 
#'   decreasing dr4pl curve as well as measured points. Default points are 
#'   blue and size 5.
#' @title plot
#' @name plot.dr4pl
#' @param x `dr4pl' object whose data and mean response function will be plotted.
#' @param text.title Character string for the title of a plot with a default set to 
#'   "Dose response plot".
#' @param text.x Character string for the x-axis of the plot with a default set to 
#'   "Dose".
#' @param text.y Character string for the y-axis of the plot with a default set to 
#'   "Response".
#' @param indices.outlier Pass a vector indicating all indices which are outliers in 
#'   the data.
#' @param ... All arguments that can normally be passed to ggplot.
#' @examples
#' ryegrass.dr4pl <- dr4pl::dr4pl(Response ~ Dose, data = sample_data_1)
#'
#' plot(ryegrass.dr4pl)
#' 
#' ##Able to further edit plots
#' library(ggplot2)
#' ryegrass.dr4pl <- dr4pl::dr4pl(Response ~ Dose, 
#'                                data = sample_data_1, 
#'                                text.title = "Sample Data Plot")
#'
#' a <- plot(ryegrass.dr4pl) 
#' a + geom_point(color = "green", size = 5)
#' 
#' ##Bring attention to outliers using parameter indices.outlier.
#' 
#' a <- dr4pl(Response ~ Dose, 
#'            data = drc_error_3, 
#'            method.init = "Mead", 
#'            method.robust = "absolute" )
#' plot(a, indices.outlier = c(90, 101))
#' 
#' ##Change the plot title default with parameter text.title
#' 
#' ryegrass.dr4pl <- dr4pl::dr4pl(Response ~ Dose, 
#'                                data = sample_data_1)
#' plot(ryegrass.dr4pl, text.title = "My New Dose Response plot")
#' 
#' ##Change the labels of the x and y axis to your need
#' 
#' library(drc)  #example requires decontaminants dataset from drc package.
#' d <- subset(decontaminants, group %in% "hpc")
#' e <- dr4pl(count~conc, data = d)
#' plot(e, 
#'      text.title = "hpc Decontaminants Plot", 
#'      text.x = "Concentration", 
#'      text.y = "Count")
#' 
#' @export
plot.dr4pl <- function(x,
                       text.title = "Dose response plot",
                       text.x = "Dose",
                       text.y = "Response",
                       indices.outlier = NULL,
                       ...) {

  ### Check whether function arguments are appropriate
  if(!is.character(text.title)) {
    
    stop("Title text should be characters.")
  }
  if(!is.character(text.x)) {
    
    stop("The x-axis label text should be characters.")
  }
  if(!is.character(text.y)) {
    
    stop("The y-axis label text should be characters.")
  }
  
  ### Draw a plot
  n <- x$sample.size
  color.vec <- rep("blue", n)
  
  if(!is.null(indices.outlier)) {
    
    color.vec[indices.outlier] <- "red"
  }
  
  a <- ggplot2::ggplot(aes(x = x$data$Dose, y = x$data$Response), data = x$data)

  a <- a + ggplot2::stat_function(fun = MeanResponse,
                                  args = list(theta = x$parameters),
                                  size = 1.2)

  a <- a + ggplot2::geom_point(size = I(5), alpha = I(0.8), color = color.vec)

  a <- a + ggplot2::labs(title = text.title,
                         x = text.x,
                         y = text.y)
  
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

  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)
}

#' @description Print the dr4pl object summary.
#' @title summary
#' @name summary.dr4pl
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

#' @description Perform the goodness-of-fit (gof) test for the 4PL model when there
#'   are at least two replicates for each dose level.
#' @title Perform the goodness-of-fit (gof) test for the 4PL model.
#' @name gof.dr4pl
#'   
#' @param object An object of the dr4pl class.
#' 
#' @return The test statistic value, its p-value and the corresponding degrees of
#'   freedom.
#'   
#' @details Perform the goodness-of-fit (gof) test for the 4PL model in which the
#'   mean response actually follws the 4 Parameter Logistic model. There should
#'   be at least two replicates at each dose level. The test statistic follows the
#'   Chi squared distributions with the (n - 4) degrees of freedom where n is the
#'   number of observations and 4 is the number of parameters in the 4PL model. For
#'   detailed explanation of the method, please refer to Subsection 2.1.5 of
#'   Seber, G. A. F. and Wild, C. J. (1989). Nonlinear Regression. Wiley Series in
#'   Probability and Mathematical Statistics: Probability and Mathematical
#'   Statistics. John Wiley & Sons, Inc., New York.
#'   
#' @author Hyowon An, Justin T. Landis and Aubrey G. Bailey
#' @export
gof.dr4pl <- function(object) {
  
  x <- object$data$Dose
  y <- object$data$Response
  
  return(NULL)
  
}

#' @description Compute the confidence intervals of parameter estimates of a fitted
#'   model.
#' @title Fit a 4 parameter logistic (4PL) model to dose-response data.
#' @name confint.dr4pl
#'   
#' @param object An object of the dr4pl class
#' @param parm Parameters of the 4PL model
#' @param level Confidence level
#' @param ... Other parameters to be passed
#' 
#' @return A matrix of the confidence intervals in which each row represents a
#'   parameter and each column represents the lower and upper bounds of the
#'   confidence intervals of the corresponding parameters.
#'   
#' @details This function computes the confidence intervals of the parameters of the
#'   4PL model based on the second order approximation to the Hessian matrix of the
#'   loss function of the model. Please refer to Subsection 5.2.2 of 
#'   Seber, G. A. F. and Wild, C. J. (1989). Nonlinear Regression. Wiley Series in
#'   Probability and Mathematical Statistics: Probability and Mathematical
#'   Statistics. John Wiley & Sons, Inc., New York.
#'   
#' @examples
#'   obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_1)
#'   parm <- obj.dr4pl$parameters
#'
#'   confint(obj.dr4pl, parm = parm, level = 0.95)
#' 
#' @author Hyowon An, Justin T. Landis and Aubrey G. Bailey
#' @export
confint.dr4pl <- function(object, parm, level, ...) {
  
  x <- object$data$Dose
  y <- object$data$Response
  theta <- object$parameters
  hessian <- object$hessian
  
  n <- object$sample.size  # Number of observations in data
  f <- MeanResponse(x, theta)
  
  C.hat.inv <- solve(hessian/2)
  s <- sqrt(sum((y - f)^2)/(n - 4))
  
  q.t <- qt(1 - (1 - level)/2, df = n - 4)
  std.err <- s*sqrt(diag(C.hat.inv))  # Standard error
  ci <- cbind(theta - q.t*std.err, theta + q.t*std.err)
  
  return(ci)
}

# add more description
#' @description Coefficient of a `dr4pl' object
#' @title coef
#' @name coef.dr4pl
#' @param object A 'dr4pl' object
#' @param ... arguments passed to coef
#' @return A vector of parameters
#' @export
coef.dr4pl <- function(object, ...) {
  
  object$parameters
  
}

