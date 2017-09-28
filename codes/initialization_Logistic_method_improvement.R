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