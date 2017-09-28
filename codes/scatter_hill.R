
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

data.drr <- read.csv("drc_error_1.csv")

drr <- drra(Response ~ Dose, data = data.drr, method.robust = "absolute")

ComparisonPlot(data.drr$Dose, data.drr$Response, coef(drr))

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Presentation")
eesetEPS()
postscript("scatter_hill_plot.eps", width = 8, height = 4)

par(mfrow = c(1, 2), mar = c(5, 4, 3, 2) + 0.1)

plot(x = log10(data.to.comp$Dose), y = data.to.comp$Response, pch = 16,
     xlab = "",
     ylab = "Response",
     main = "Scatter plot",
     bty = "L",
     xaxt = "n")

mtext(text = "Dose", side = 1, line = 4)

x.levels <- sort(unique(data.to.comp$Dose))
at.x <- log10(x.levels)
axis(1, at = at.x, labels = 10^at.x, las = 2)



dev.off()
