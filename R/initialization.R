
#' Find initial values for the 4PL model.
#'
#' Args:
#' @param x x values
#' @param y y values
#' @param method.init Initialization method
#' @param method.robust Robust fitting method
#'
#' @returns theta.IC50.slope Parameter estimates of the IC50 and slope
#' @export
FindInitialParms <- function(x, y, method.init, method.robust) {

  scale.inc <- 0.001
  y.range <- range(y)
  len.y.range <- scale.inc * diff(y.range)

  theta.1.init <- max(y) + len.y.range
  theta.4.init <- min(y) - len.y.range

  # If theta.1 < theta.4, then the curve is a growth curve.
  # If theta.1 > theta.4, then the curve is a decline curve.

  if(method.init == "logistic") {

    y.transf <- log10((y - theta.4.init)/(theta.1.init - y))
    x.log10 <- log10(x)

    data.lm <- data.frame(x = x.log10, y = y.transf)
    data.lm <- data.lm[data.lm$x != -Inf, ]

    #lm.init <- stats::lm(y ~ x, data = data.lm)  # Linear model for initial parameter estimates
    lm.init <- stats::lm(y ~ x, data = data.lm)  # Linear model for initial parameter estimates
    beta.hat <- lm.init$coefficients

    theta.3.init <- beta.hat[2]
    theta.2.init <- 10^(-beta.hat[1]/theta.3.init)

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

