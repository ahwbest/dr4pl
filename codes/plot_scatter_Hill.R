

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

file.name <- "Dirk_data_9.csv"
df.drc <- read.csv(file.name)
df.drc <- df.drc[df.drc$Dose != 0, ]

#df.drc <- read.csv(".csv")
#df.drc <- df.drc[df.drc$Dose != 0, ]

x.orig.vec <- df.drc$Dose
y.orig.vec <- df.drc$Response

### Scatter plot and contour plot
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Presentation")
setEPS()
postscript("hill_plot_Dirk_data_9.eps", width = 8, height = 4)

par(mfrow = c(1, 2), mar = c(5, 4, 3, 2) + 0.1)

# Scatter plot
plot(x = log10(x.orig.vec),
     y = y.orig.vec,
     main = "Scatter plot",
     xlab = "",
     ylab = "Response",
     lwd = 2,
     pch = 16,
     type = "p",
     bty = "L",
     xaxt = "n")

mtext(text = "Dose", side = 1, line = 4)

x.levels <- sort(unique(x.orig.vec))
at.x <- log10(x.levels)
axis(1, at = at.x, labels = 10^at.x, las = 2)

# Hill plot
x <- x.orig.vec
y <- y.orig.vec

scale.inc <- 0.1
y.range <- range(y)
len.y.range <- scale.inc * diff(y.range)

y.max <- max(y) + len.y.range
y.min <- min(y) - len.y.range

y.logit <- log((y - y.min)/(y.max - y))  # Logit transformed responses
x.log <- log(x)  # Log transformed doses

lm.hill <- lm(y.logit ~ x.log)
rlm.hill <- ltsReg(y.logit ~ x.log)
beta.lm <- coef(lm.hill)
beta.rlm <- coef(rlm.hill)

plot(x = x.log, y = y.logit, pch = 16)
abline(a = beta.lm[1], b = beta.lm[2], lwd = 2, col = "blue")
abline(a = beta.rlm[1], b = beta.rlm[2], lwd = 2, col = "red")
legend(x = "topright", col = c("blue", "red"), lwd = 2,
       legend = c("Nonrobust", "Robust"))

dev.off()

# theta.3.init <- beta.lm[2]
# theta.2.init <- exp(-beta.lm[1]/theta.3.init)
# theta.1.init <- y.max
# theta.4.init <- y.min
# 
# theta.init <- c(theta.1.init, theta.2.init, theta.3.init, theta.4.init)
# 
# obj.drc <- drm(Response ~ Dose,
#                data = df.drc,
#                fct = LL.4(names = c("Slope", "Left limit", "Right limit", "IC50")),
#                start = theta.init)

