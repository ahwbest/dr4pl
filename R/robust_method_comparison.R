# -----------------------------------------------------------------------------
### Load libraries and source code
#
library(drc)
source(".\\R\\initialization.R")
source(".\\R\\main.R")

# -----------------------------------------------------------------------------
### Test functions for initialization
#
ComparisonPlot <- function(dose, response,
                           parm.mat,
                           legend.names = NULL) {
  # Make a scatter plot with dose response curves drawn on it.
  #
  # Args:
  #   dose, response: x and y values to be plotted on the scatter plot.
  #   parm.mat: The matrix of parameters. Each row represents the parameters of
  #             a single fitted model.
  #   legend.names: The names of legends. If NULL, no legend is displayed.
  #
  # Returns:
  #   Display a plot. No return value.
  data.drm <- data.frame(Dose = dose, Response = response)

  a <- ggplot(aes(x = Dose, y = Response), data = data.drm)

  for(i in 1:nrow(parm.mat)) {
    a <- a + stat_function(fun = MeanResponseCurve,
                           args = list(theta = parm.mat[i, ]),
                           colour = i,
                           size = 1.2)
  }

  a <- a + geom_point(size = I(5), alpha = I(0.8))

  a <- a + labs(title = "Dose response curve")

  # Set parameters for the grids
  a <- a + theme(strip.text.x = element_text(size = 16))
  a <- a + theme(panel.grid.minor = element_blank())
  a <- a + theme(panel.grid.major = element_blank())
  a <- a + scale_x_log10()
  a <- a + theme_bw()

  # Set parameters for the titles and text / margin(top, right, bottom, left)
  a <- a + theme(plot.title = element_text(size = 20, margin = margin(0, 0, 10, 0)))
  a <- a + theme(axis.title.x = element_text(size = 16, margin = margin(15, 0, 0, 0)))
  a <- a + theme(axis.title.y = element_text(size = 16, margin = margin(0, 15, 0, 0)))
  a <- a + theme(axis.text.x = element_text(size = 16))
  a <- a + theme(axis.text.y = element_text(size = 16))

  # if(!is.null(lab.names)) {
  #   if(length(lab.names) != nrow(parm.mat)) {
  #     stop("The number of labes is not equal to the number of curves.")
  #   }
  #
  #   a <- a + theme(legend.title = element_text(size = 14))
  #   a <- a + theme(legend.text = element_text(size = 12))
  #   a <- a + theme(legend.key = element_rect(colour = "white"))
  # }

  plot(a)
}

CompareRobustMethods <- function(data.to.comp,
                                 ind.plot = FALSE,
                                 var.dose,
                                 var.response,
                                 var.ref = NULL) {
  # Compare different initialization methods.
  #
  # Args:
  #   data.to.comp: data to fit both drm and drra
  #   var.dose: variable corresponding to dose
  #   var.response: variable corresponding to response
  #   var.ref: reference variable
  #
  # Returns:
  #   Print the estimated parameters to standard output
  data.to.comp <- na.omit(data.to.comp)

  if(length(var.ref) == 0) {

    data.new <- subset(x = data.to.comp, select = c(var.dose, var.response))
    colnames(data.new) <- c("Dose", "Response")
    n <- nrow(data.new)  # Number of observations

    # Four different robust methods
    obj.drra.1 <- drra(Response ~ Dose, data = data.new)
    obj.drra.2 <- drra(Response ~ Dose, data = data.new,
                       method.robust = "absolute")
    obj.drra.3 <- drra(Response ~ Dose, data = data.new,
                       method.robust = "Huber")
    obj.drra.4 <- drra(Response ~ Dose, data = data.new,
                       method.robust = "Tukey")

    parm.drra.1 <- coef(obj.drra.1)
    parm.drra.2 <- coef(obj.drra.2)
    parm.drra.3 <- coef(obj.drra.3)
    parm.drra.4 <- coef(obj.drra.4)

    if(ind.plot) {
      x <- data.new$Dose
      y <- data.new$Response

      plot(x = x, y = y
           pch = 16)
    }
  }
}

# -------------------------------------------------------------------------------
### Compare the parameter estimates
#
### algae
cat("algae\n")

# Plot fitted curves
CompareRobustMethods(data.to.comp = algae,
                     ind.plot = TRUE,
                     var.dose = "conc",
                     var.response = "vol")

### etmotc
cat("etmotc\n")

CompareInitializationMethods(data.to.comp = etmotc,
                             var.ref = "pct1",
                             var.dose = "dose1",
                             var.response = "rgr1")

### G.aparine
cat("G.aparine\n")

CompareInitializationMethods(data.to.comp = G.aparine,
                             var.dose = "dose",
                             var.response = "drymatter")

### glymet
cat("glylmet\n")

CompareInitializationMethods(data.to.comp = glymet,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr")

### heartrate
cat("heartrate\n")

CompareInitializationMethods(data.to.comp = heartrate,
                             var.dose = "pressure",
                             var.response = "rate")

### leaflength
cat("leaflength\n")

CompareInitializationMethods(data.to.comp = leaflength,
                             var.dose = "Dose",
                             var.response = "DW")

### lepidium
cat("lepidium\n")

CompareInitializationMethods(data.to.comp = lepidium,
                             var.dose = "conc",
                             var.response = "weight")

### lettuce
cat("lettuce\n")

CompareInitializationMethods(data.to.comp = lettuce,
                             var.dose = "conc",
                             var.response = "weight")

### mecter
cat("mecter\n")

CompareInitializationMethods(data.to.comp = mecter,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr")

### M.bahia
cat("M.bahia\n")

CompareInitializationMethods(data.to.comp = M.bahia,
                             var.dose = "conc",
                             var.response = "dryweight")

### nasturtium
cat("nasturtium\n")

CompareInitializationMethods(data.to.comp = nasturtium,
                             var.dose = "conc",
                             var.response = "weight")

### O.mykiss
cat("O.mykiss\n")

CompareInitializationMethods(data.to.comp = O.mykiss,
                             var.dose = "conc",
                             var.response = "weight")

### P.promelas
cat("P.promelas\n")

CompareInitializationMethods(data.to.comp = P.promelas,
                             var.dose = "conc",
                             var.response = "dryweight")

### ryegrass
cat("ryegrass\n")

CompareInitializationMethods(data.to.comp = ryegrass,
                             var.dose = "conc",
                             var.response = "rootl")

### S.alba
cat("S.alba\n")

CompareInitializationMethods(data.to.comp = S.alba,
                             var.ref = "Herbicide",
                             var.dose = "Dose",
                             var.response = "DryMatter")

### S.capricornutum
cat("S.capricornutum\n")

CompareInitializationMethods(data.to.comp = S.capricornutum,
                             var.dose = "conc",
                             var.response = "count")

### secalonic
cat("secalonic\n")

CompareInitializationMethods(data.to.comp = secalonic,
                             var.dose = "dose",
                             var.response = "rootl")

### spinach
cat("spinach\n")

CompareInitializationMethods(data.to.comp = spinach,
                             var.ref = "CURVE",
                             var.dose = "DOSE",
                             var.response = "SLOPE")

### terbuthylazin
cat("terbuthylazin\n")

CompareInitializationMethods(data.to.comp = terbuthylazin,
                             var.dose = "dose",
                             var.response = "rgr")

### vinclozolin
cat("vinclozolin\n")

CompareInitializationMethods(data.to.comp = vinclozolin,
                             var.ref = "exper",
                             var.dose = "conc",
                             var.response = "effect")

sink()
