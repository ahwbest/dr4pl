
#' Find initial parameter estimates for the left and right asymptotes in the
#' 4PL model.
#'
#' @param x x values
#' @param y y values
#'
#' @return theta.left.right: Parameter estimates of the left and right asymptotes.
#' @export
FindLeftRightAsymptotes <- function(x, y) {

  scale.inc <- 0.001
  y.range <- range(y)
  len.y.range <- scale.inc * diff(y.range)

  theta.1.init <- max(y) + len.y.range
  theta.4.init <- min(y) - len.y.range

  return(c(theta.1.init, theta.4.init))
}


#' Find initial values for the IC50 and slope parameters.
#' @param x x values
#' @param y y values
#' @param theta.1.4.init Estimates of the left and right asymptotes
#' @param method.init Initialization method
#'
#' @return theta.IC50.slope: Parameter estimates of the IC50 and slope
#' @export
FindIC50Slope <- function(x, y, theta.1.4.init, method.init) {

  theta.1.init <- theta.1.4.init[1]
  theta.4.init <- theta.1.4.init[2]

  if(method.init == "logistic") {

    y.transf <- log10((y - theta.4.init)/(theta.1.init - y))
    x.log10 <- log10(x)

    data.lm <- data.frame(x = x.log10, y = y.transf)
    data.lm <- data.lm[data.lm$x != -Inf, ]

    lm.init <- stats::lm(y ~ x, data = data.lm)  # Linear model for initial parameter estimates
    beta.hat <- lm.init$coefficients

    theta.3.init <- beta.hat[2]
    theta.2.init <- 10^(-beta.hat[1]/theta.3.init)

  } else if(method.init == "Mead") {

    log.x <- log(x)
    y.zero.low <- y - theta.4.init

    gamma.seq <- seq(from = 0.05, to = 0.95, by = 0.05)

    for(gamma in gamma.seq) {
      mat <- matrix(nrow = 2, ncol = 2)
      mat[1, 1] <- sum(y.zero.low^2)
      mat[2, 1] <- mat[1, 2] <- gamma^log.x%*%y.zero.low^2
      mat[2, 2] <- gamma^(2*log.x)%*%y.zero.low^2

      vec <- c(sum(y.zero.low), gamma^log.x%*%y.zero.low)
    }
  }

  return(c(theta.2.init, theta.3.init))
}

