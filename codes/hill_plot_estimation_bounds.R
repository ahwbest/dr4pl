setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Presentation")
setEPS()
postscript("hill_plot_estimation_bounds.eps", width = 8, height = 4)

par(mfrow = c(1, 2), mar = c(5, 4, 3, 2) + 0.1)

# Hill plot
scale.inc <- 0.001
y.range <- range(y)
len.y.range <- scale.inc * diff(y.range)

y.max <- max(y) + len.y.range
y.min <- min(y) - len.y.range

y.logit <- log((y - y.min)/(y.max - y))  # Logit transformed responses
x.log <- log(x)  # Log transformed doses

lm.hill <- lm(y.logit ~ x.log)
beta.lm <- coef(lm.hill)

plot(x = x.log, y = y.logit, pch = 16)
abline(a = beta.lm[1], b = beta.lm[2], lwd = 2)

# Hill plot
scale.inc <- 0.001
y.range <- range(y)
len.y.range <- scale.inc * diff(y.range)

y.max <- max(y) + len.y.range
y.min <- min(y) - len.y.range

y.logit <- log((y - theta.4)/(y.max - y))  # Logit transformed responses
x.log <- log(x)  # Log transformed doses

lm.hill <- lm(y.logit ~ x.log)
beta.lm <- coef(lm.hill)

plot(x = x.log, y = y.logit, pch = 16)
abline(a = beta.lm[1], b = beta.lm[2], lwd = 2)

dev.off()