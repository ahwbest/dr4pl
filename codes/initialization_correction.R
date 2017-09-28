# -----------------------------------------------------------------------------
### Load libraries and source codes
#
library(drc)
source(".\\R\\initialization.R")
source(".\\R\\main.R")

# -----------------------------------------------------------------------------
### User defined functions
#
findbe1 <- function(doseTr, respTr, sgnb = 1, back = exp)
{
  function(x, y, cVal, dVal)
  {
    #        lmFit <- lm(log((dVal - y)/(y - cVal)) ~ log(x), subset = x > 0)

    #respTr <- function(y, cVal, dVal) {log((dVal - y)/(y - cVal))}
    #doseTr <- function(x) {rVec <- log(x); rVec[!x>0] <- NA; rVec}
    lmFit <- lm(respTr(y, cVal, dVal) ~ doseTr(x))
    coefVec <- coef(lmFit)
    bVal <- sgnb * coefVec[2]
    eVal <- back(-coefVec[1] / (sgnb * bVal))

    return(as.vector(c(bVal, eVal)))
  }
}

findcd <- function(x, y, scaleInc = 0.001)
{
  yRange <- range(y)
  lenyRange <- scaleInc * diff(yRange)
  #    cVal <- yRange[1] - lenyRange  # the c parameter
  #    dVal <- yRange[2] + lenyRange  # the d parameter

  c(yRange[1] - lenyRange, yRange[2] + lenyRange)
}

FindLeftRightAsymptotes <- function(x, y) {
  # Find initial parameter estimates for the left and right asymptotes in the
  # 4PL model.
  #
  # Args:
  #   x: x values
  #   y: y values
  #
  # Returns:
  #   theta.left.right: Parameter estimates of the left and right asymptotes.
  scale.inc <- 0.001
  y.range <- range(y)
  len.y.range <- scale.inc * diff(y.range)

  theta.1.init <- max(y) + len.y.range
  theta.4.init <- min(y) - len.y.range

  return(c(theta.1.init, theta.4.init))
}

FindIC50Slope <- function(x, y, theta.1.4.init, method.init) {
  # Find initial values for the IC50 and slope parameters.
  #
  # Args:
  #   x: x values
  #   y: y values
  #   theta.1.4.init: Estimates of the left and right asymptotes
  #   method.init: Initialization method
  #
  # Returns:
  #   theta.IC50.slope: Parameter estimates of the IC50 and slope
  theta.1.init <- theta.1.4.init[1]
  theta.4.init <- theta.1.4.init[2]

  if(method.init == "logistic") {
    y.transf <- log10((y - theta.4.init)/(theta.1.init - y))
    x.log10 <- log10(x)

    data.lm <- data.frame(x = x.log10, y = y.transf)
    data.lm <- data.lm[data.lm$x != -Inf, ]

    lm.init <- lm(y ~ x, data = data.lm)  # Linear model for initial parameter estimates
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

"llogistic.ssf" <- function(method = c("1", "2", "3", "4"), fixed, useFixed = FALSE)
{
  method <- match.arg(method)

  ## Defining helper functions (used below)
  ytrans <- function(y, cVal, dVal) {log((dVal - y)/(y - cVal))}
  bfct <- function(x, y, cVal, dVal, eVal) {ytrans(y, cVal, dVal) / log(x / eVal)}
  #    efct <- function(x, y, bVal, cVal, dVal) {x * (((dVal - y) / (y - cVal))^(-1 / bVal))}
  efct <- function(x, y, bVal, cVal, dVal) {x * exp(-ytrans(y, cVal, dVal)/bVal)}

  ## Assigning function for finding initial b and e parameter values
  findbe <- switch(method,
                   "1" = findbe1(function(x) {rVec <- log(x); rVec[!x>0 | !is.finite(x)] <- NA; rVec}, ytrans),
                   "2" = findbe2(bfct, efct, "Anke"),
                   "3" = findbe3(),
                   "4" = findbe2(bfct, efct, "Normolle"))

  function(dframe)
  {
    ncoldf <- ncol(dframe)
    x <- dframe[, 1]
    #        x <- dframe[, -ncoldf]
    y <- dframe[, ncoldf]

    #        x <- dframe[, 1]
    #        y <- dframe[, 2]

    ## Finding initial values for c and d parameters
    cdVal <- findcd(x, y)
    #        if (useFixed) {  # not implemented at the moment
    #            cdVal <- c(ifelse(notFixed[2], cdVal[1], fixed[2] / respScaling),
    #            ifelse(notFixed[3], cdVal[2], fixed[3] / respScaling))}

    ## Finding initial values for b and e parameters
    beVal <- findbe(x, y, cdVal[1], cdVal[2])

    ## Finding initial value for f parameter
    fVal <- 1
    # better choice than 1 may be possible!
    # the f parameter, however, is very rarely a magnitude of 10 larger or smaller

    return(c(beVal[1], cdVal, beVal[2], fVal)[is.na(fixed)])
  }
}

# -----------------------------------------------------------------------------
### Main part
#
# data.orgn = glymet
# var.ref = "pct"
# var.dose = "dose"
# var.response = "rgr"
#
# data.target <- subset(x = data.orgn, select = c(var.ref, var.dose, var.response))
# colnames(data.target) <- c("Ref", "Dose", "Response")
# data.target$Ref <- as.factor(data.target$Ref)
# data.target <- subset(x = data.target, subset = Ref != "999")
# data.target <- droplevels(data.target)
#
# levels.ref <- levels(data.target$Ref)
# ref <- "67"
#
# data.drm <- subset(x = data.target, select = c(Dose, Response), subset = Ref == ref)
# n <- nrow(data.drm)  # Number of observations
# x <- data.drm$Dose
# y <- data.drm$Response
#
# # Obtain the starting initial values of drc
# ssf.drc <- llogistic.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA), useFixed = FALSE)
# parm.init.drc <- ssf.drc(data.drm)
# parm.init.drc <- parm.init.drc[-5]
#
# # Obtain the starting initial values of drra
# theta.1.4.init <- FindLeftRightAsymptotes(x, y)
# theta.2.3.init <- FindIC50Slope(x, y, theta.1.4.init, method.init = "logistic")
# parm.init.drra <- c(theta.2.3.init[2], theta.1.4.init[1], theta.1.4.init[2], theta.2.3.init[1])
#
# obj.drc.1 <- drm(Response ~ Dose,
#                  data = data.drm,
#                  fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
#                             method = "1"),
#                  control = drmc(method = "Nelder-Mead"),
#                  start = parm.init.drc)
#
# obj.drc.2 <- drm(Response ~ Dose,
#                  data = data.drm,
#                  fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
#                             method = "1"),
#                  control = drmc(method = "Nelder-Mead"),
#                  start = parm.init.drra)
#
# obj.drc.3 <- drm(Response ~ Dose,
#                  data = data.drm,
#                  fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
#                             method = "1"),
#                  control = drmc(method = "Nelder-Mead"))
#
# parm.drc <- rbind(coef(obj.drc.1), coef(obj.drc.2), coef(obj.drc.3))
#
# obj.drra <- drra(Response ~ Dose,
#                  data = data.each)
#
# parm.drra <- coef(obj.drra)
