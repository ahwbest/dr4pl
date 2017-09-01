
#' Find initial values for the 4PL model.
#'
#' Args:
#' @param x x values
#' @param y y values
#' @param method.init Initialization method
#' @param method.robust Robust fitting method
#'
#' @return theta.IC50.slope Parameter estimates of the IC50 and slope
#' @export
FindInitialParms <- function(x, y, method.init, method.robust) {

  scale.inc <- 1  #Should this value be 0.001 according to what is on your paper??
  y.range <- range(y)
  len.y.range <- scale.inc * diff(y.range)

  y.max <- max(y) + len.y.range
  y.min <- min(y) - len.y.range

  ### Check whether input are appropriate.
  if(length(x) == 0 || length(y) == 0 || length(x) != length(y)) {

    stop("The same numbers of dose and response values should be supplied.")
  }

  ### If theta.1 < theta.4, then the curve is a growth curve.
  ### If theta.1 > theta.4, then the curve is a decline curve.
  if(method.init == "logistic") {

    y.logit <- log((y - y.min)/(y.max - y))  # Logit transformed responses
    x.log <- log(x)  # Log transformed doses

    data.hill <- data.frame(y.logit = y.logit,
                            x.log = x.log,
                            y.max.y = 1/(y.max - y),
                            y.y.min = 1/(y - y.min))

    data.hill <- subset(data.hill, subset = x.log != -Inf) #should take out any x=0

    lm.hill <- lm(y.logit ~ x.log + y.max.y + y.y.min, data = data.hill) #paper makes no mention of y.max.y and y.y.min in the regression model
    lm.hill.simple <- lm(y.logit ~ x.log, data = data.hill)

    # summary(lm.hill)
    # summary(lm.hill.simple)

    # plot(x = data.hill$x.log, y = data.hill$y.logit, pch = 16)
    # abline(a = lm.hill$coefficients[1], b = lm.hill$coefficients[2],
    #        col = "blue", lwd = 2)
    # plot(x = data.hill$x.log, y = data.hill$y.logit - data.hill)
    # abline(a = lm.hill.simple$coefficients[1], b = lm.hill.simple$coefficients[2],
    #        col = "red", lwd = 2)

    # Fully estimate the IC50 and slope parameters first, and then others
    data.resid <- data.frame(res = lm.hill.simple$residuals,
                             y.max.y = data.hill$y.max.y,
                             y.y.min = data.hill$y.y.min)
    lm.resid <- lm(res ~ y.max.y + y.y.min, data = data.resid)

    # plot(x = log10(x), y = y, pch = 16)

    #lm.hill <- ltsReg(y.logit ~ x.log + y.max.y + y.y.min, data = data.hill)
    #lm.hill <- lmrob(y.logit ~ x.log + y.max.y + y.y.min,
    #                 data = data.hill,
    #                 maxit.scale = 500)

    # beta.hat <- lm.hill$coefficients
    #
    # theta.1.init <- beta.hat[3]
    # theta.3.init <- beta.hat[2]
    # theta.2.init <- exp(-beta.hat[1]/theta.3.init)
    # theta.4.init <- beta.hat[4]

    beta.hat <- lm.hill.simple$coefficients
    theta.1.init <- y.max
    theta.3.init <- beta.hat[2]
    theta.2.init <- exp(-beta.hat[1]/theta.3.init)
    theta.4.init <- y.min

    theta.init <- c(theta.1.init, theta.2.init, theta.3.init, theta.4.init)
  } else if(method.init == "Mead") {

    log.x <- log(x)
    y.lower.bd <- min(theta.1.init, theta.4.init)
    y.zero.low <- y - y.lower.bd

    gammas.temp <- seq(from = 0.05, to = 0.95, by = 0.05)  # Reciprocals of gamma values
    gammas <- c(gammas.temp, rev(1/gammas.temp))
    theta.matr <- matrix(0, nrow = length(gammas), ncol = 4)

    for(i in 1:length(gammas)) {

      gamma <- gammas[i]

      data.lm <- data.frame(y = y.zero.low,
                            y.gamma.x = y.zero.low*gamma^log.x,
                            Response = 1)

      # Set the second predictor values to be zero when dose is zero
      data.lm$y.gamma.x[x == 0] <- 0

      lm.init <- stats::lm(Response ~ -1 + y + y.gamma.x, data = data.lm)

      alpha.beta <- lm.init$coefficients
      alpha <- alpha.beta[1]
      beta <- alpha.beta[2]

      theta.matr[i, 1] <- 1/alpha

      if(alpha*beta > 0) {

        theta.matr[i, 2] <- exp((log(alpha/beta))/log(gamma))
      } else {

        theta.matr[i, 2] <- 0
      }

      theta.matr[i, 3] <- log(gamma)
      theta.matr[i, 4] <- 0
    }

    theta.matr[, 1] <- theta.matr[, 1] + y.lower.bd
    theta.matr[, 4] <- theta.matr[, 4] + y.lower.bd

    err.fcn <- ErrFcn(method.robust)
    errors <- rep(0, length(gammas))  # To save the error values

    for(i in 1:length(gammas)) {

      theta <- theta.matr[i, ]
      errors[i] <- err.fcn(theta, x, y)
    }

    theta.init <- theta.matr[which.min(errors), ]
  }

  names(theta.init) <- c("theta.1", "theta.2", "theta.3", "theta.4")

  return(theta.init)
}

