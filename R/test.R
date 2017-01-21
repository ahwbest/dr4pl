# ------------------------------------------------------------------------
### Load libraries
#
library(drc)

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

# ------------------------------------------------------------------------
### Error case 1
#
data.error.1 <- read.csv("drc_error_1.csv")

# ------------------------------------------------------------------------
### Generate error case 2
#
x <- rep(c(0.00135, 0.0135, 0.135, 1.35, 13.5), each = 2)
y <- c(10, 20, 18, 22, 30, 70, 10, 100, 5, 7)

data.error.2 <- data.frame(Dose = x, Response = y)

# ------------------------------------------------------------------------
###
#
data.error.3 <- read.csv("Dirk_data_inh4_80.csv")

data.error.3 <- subset(subset = Plate == 37,
                       select = c(RLU, Conc1),
                       x = data.error.3)

colnames(data.error.3) <- c("Response", "Dose")

drmc.error <- drmc(method = "Nelder-Mead")

drm.error.1 <- drm(Response ~ Dose,
                   control = drmc.error,
                   data = data.error.3,
                   fct = LL.4())

drm.error.2 <- drm(Response ~ Dose,
#                   control = drmc.error,
                   data = data.error.2,
                   fct = LL.4())

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Submission_JSS")
jpeg(filename = "drc_error_cases.jpg",
     width = 8, height = 4, units = "in", res = 600)

par(mfrow = c(1, 2), mar = c(5, 4, 3, 1) + 0.1)
plot(x = log10(data.error.1$Dose), y = data.error.1$Response,
     pch = 16,
     main = "Error Case 1",
     xlab = expression("log"[10]*"(dose)"),
     ylab = "Response",
     bty = "l",
     xaxt = "n")

axis(side = 1, at = log10(unique(data.error.1$Dose)),
     labels = unique(data.error.1$Dose))

plot(x = log10(data.error.2$Dose), y = data.error.2$Response,
     pch = 16,
     main = "Error Case 2",
     xlab = expression("log"[10]*"(dose)"),
     ylab = "Response",
     bty = "l",
     xaxt = "n")

axis(side = 1, at = log10(unique(data.error.2$Dose)),
     labels = unique(data.error.2$Dose))

dev.off()
