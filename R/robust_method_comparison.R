# -----------------------------------------------------------------------------
### Load libraries and source code
#
library(drc)
source(".\\R\\initialization.R")
source(".\\R\\main.R")

# -----------------------------------------------------------------------------
### Test functions for initialization
#
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

      plot(x = x, y = y, pch = 16)
    }
  }
}

# -------------------------------------------------------------------------------
### Compare the parameter estimates
#
### algae
# cat("algae\n")
#
# # Plot fitted curves
# CompareRobustMethods(data.to.comp = algae,
#                      ind.plot = TRUE,
#                      var.dose = "conc",
#                      var.response = "vol")
#
# ### etmotc
# cat("etmotc\n")
#
# CompareRobustMethods(data.to.comp = etmotc,
#                              var.ref = "pct1",
#                              var.dose = "dose1",
#                              var.response = "rgr1")
#
# ### G.aparine
# cat("G.aparine\n")
#
# CompareRobustMethods(data.to.comp = G.aparine,
#                              var.dose = "dose",
#                              var.response = "drymatter")
#
# ### glymet
# cat("glylmet\n")
#
# CompareRobustMethods(data.to.comp = glymet,
#                              var.ref = "pct",
#                              var.dose = "dose",
#                              var.response = "rgr")
#
# ### heartrate
# cat("heartrate\n")
#
# CompareRobustMethods(data.to.comp = heartrate,
#                              var.dose = "pressure",
#                              var.response = "rate")
#
# ### leaflength
# cat("leaflength\n")
#
# CompareRobustMethods(data.to.comp = leaflength,
#                              var.dose = "Dose",
#                              var.response = "DW")
#
# ### lepidium
# cat("lepidium\n")
#
# CompareRobustMethods(data.to.comp = lepidium,
#                              var.dose = "conc",
#                              var.response = "weight")
#
# ### lettuce
# cat("lettuce\n")
#
# CompareRobustMethods(data.to.comp = lettuce,
#                              var.dose = "conc",
#                              var.response = "weight")
#
# ### mecter
# cat("mecter\n")
#
# CompareRobustMethods(data.to.comp = mecter,
#                              var.ref = "pct",
#                              var.dose = "dose",
#                              var.response = "rgr")
#
# ### M.bahia
# cat("M.bahia\n")
#
# CompareRobustMethods(data.to.comp = M.bahia,
#                              var.dose = "conc",
#                              var.response = "dryweight")
#
# ### nasturtium
# cat("nasturtium\n")
#
# CompareRobustMethods(data.to.comp = nasturtium,
#                              var.dose = "conc",
#                              var.response = "weight")
#
# ### O.mykiss
# cat("O.mykiss\n")
#
# CompareRobustMethods(data.to.comp = O.mykiss,
#                              var.dose = "conc",
#                              var.response = "weight")
#
# ### P.promelas
# cat("P.promelas\n")
#
# CompareRobustMethods(data.to.comp = P.promelas,
#                              var.dose = "conc",
#                              var.response = "dryweight")
#
# ### ryegrass
# cat("ryegrass\n")
#
# CompareRobustMethods(data.to.comp = ryegrass,
#                              var.dose = "conc",
#                              var.response = "rootl")
#
# ### S.alba
# cat("S.alba\n")
#
# CompareRobustMethods(data.to.comp = S.alba,
#                              var.ref = "Herbicide",
#                              var.dose = "Dose",
#                              var.response = "DryMatter")
#
# ### S.capricornutum
# cat("S.capricornutum\n")
#
# CompareRobustMethods(data.to.comp = S.capricornutum,
#                              var.dose = "conc",
#                              var.response = "count")
#
# ### secalonic
# cat("secalonic\n")
#
# CompareRobustMethods(data.to.comp = secalonic,
#                              var.dose = "dose",
#                              var.response = "rootl")
#
# ### spinach
# cat("spinach\n")
#
# CompareRobustMethods(data.to.comp = spinach,
#                              var.ref = "CURVE",
#                              var.dose = "DOSE",
#                              var.response = "SLOPE")
#
# ### terbuthylazin
# cat("terbuthylazin\n")
#
# CompareRobustMethods(data.to.comp = terbuthylazin,
#                              var.dose = "dose",
#                              var.response = "rgr")
#
# ### vinclozolin
# cat("vinclozolin\n")
#
# CompareRobustMethods(data.to.comp = vinclozolin,
#                              var.ref = "exper",
#                              var.dose = "conc",
#                              var.response = "effect")
#
# sink()
