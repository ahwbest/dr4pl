# -----------------------------------------------------------------------------
### Load libraries and source code
#
library(drc)
library(testthat)

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

  results <- matrix(nrow = 4, ncol = 5)

  # There is only one drug
  if(length(var.ref) == 0) {

    data.each <- subset(x = data.to.comp, select = c(var.dose, var.response))
    colnames(data.each) <- c("Dose", "Response")
    n <- nrow(data.each)  # Number of observations

    # Various methods
    obj.drra <- drra(Response ~ Dose,
                     data = data.each)

    results[1, ] <- c(coef(obj.drra), obj.drra$error.value)

    obj.drra <- drra(Response ~ Dose,
                     data = data.each,
                     method.robust = "absolute")

    results[2, ] <- c(coef(obj.drra), obj.drra$error.value)

    obj.drra <- drra(Response ~ Dose,
                     data = data.each,
                     method.robust = "Huber")

    results[3, ] <- c(coef(obj.drra), obj.drra$error.value)

    obj.drra <- drra(Response ~ Dose,
                     data = data.each,
                     method.robust = "Tukey")

    results[4, ] <- c(coef(obj.drra), obj.drra$error.value)

    row.names(results) <- c("Squared", "Absolute", "Huber", "Tukey")
    colnames(results) <- c("Lower limit", "IC50", "Slope", "Upper limit", "Loss value")

    print(results)
    cat("\n")

    # Temporarily disabled
    # if(ind.plot) {
    #   # Make a plot of fitted curves
    #   parm.mat <- rbind(coef(obj.drc.1), coef(obj.drc.2), coef(obj.drc.3))
    #   parm.mat <- parm.mat[, c(3, 4, 1, 2)]
    #   parm.mat[, 3] <- -parm.mat[, 3]
    #   parm.mat <- rbind(parm.mat, parm.drra)
    #
    #   row.names(parm.mat) <- c()
    #   colnames(parm.mat) <- c("Left limit", "IC50", "Slope", "Right limit")
    #
    #   ComparisonPlot(dose = data.each$Dose, response = data.each$Response, parm.mat)
    # }

  } else {  # There are multiple drugs

    data.whole <- subset(x = data.to.comp, select = c(var.ref, var.dose, var.response))
    colnames(data.whole) <- c("Ref", "Dose", "Response")
    data.whole$Ref <- as.factor(data.whole$Ref)
    data.whole <- subset(x = data.whole, subset = Ref != "999")
    data.whole <- droplevels(data.whole)

    levels.ref <- levels(data.whole$Ref)

    for(i in 1:length(levels.ref)) {

      data.each <- subset(x = data.whole, subset = Ref == levels.ref[i])
      n <- nrow(data.each)  # Number of observations

      cat(var.ref, "=", levels.ref[i], "\n")

      # Various methods
      obj.drra <- drra(Response ~ Dose,
                       data = data.each)

      results[1, ] <- c(coef(obj.drra), obj.drra$error.value)

      obj.drra <- drra(Response ~ Dose,
                       data = data.each,
                       method.robust = "absolute")

      results[2, ] <- c(coef(obj.drra), obj.drra$error.value)

      obj.drra <- drra(Response ~ Dose,
                       data = data.each,
                       method.robust = "Huber")

      results[3, ] <- c(coef(obj.drra), obj.drra$error.value)

      obj.drra <- drra(Response ~ Dose,
                       data = data.each,
                       method.robust = "Tukey")

      results[4, ] <- c(coef(obj.drra), obj.drra$error.value)

      row.names(results) <- c("Squared", "Absolute", "Huber", "Tukey")
      colnames(results) <- c("Lower limit", "IC50", "Slope", "Upper limit", "Loss value")

      print(results)
      cat("\n")
    }
  }
}

# -------------------------------------------------------------------------------
### Compare the parameter estimates
#
sink(file = "robust_estimation_comparison.Rout")

### acidiq
cat("acidiq\n")

CompareRobustMethods(data.to.comp = acidiq,
                     var.ref = "pct",
                     var.dose = "dose",
                     var.response = "rgr")

### algae
cat("algae\n")

# Plot fitted curves
CompareRobustMethods(data.to.comp = algae,
                     var.dose = "conc",
                     var.response = "vol")

### etmotc
cat("etmotc\n")

CompareRobustMethods(data.to.comp = etmotc,
                             var.ref = "pct1",
                             var.dose = "dose1",
                             var.response = "rgr1")

### G.aparine
cat("G.aparine\n")

CompareRobustMethods(data.to.comp = G.aparine,
                             var.ref = "treatment",
                             var.dose = "dose",
                             var.response = "drymatter")

### glymet
cat("glylmet\n")

CompareRobustMethods(data.to.comp = glymet,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr")

### heartrate
cat("heartrate\n")

CompareRobustMethods(data.to.comp = heartrate,
                             var.dose = "pressure",
                             var.response = "rate")

### leaflength
cat("leaflength\n")

CompareRobustMethods(data.to.comp = leaflength,
                             var.dose = "Dose",
                             var.response = "DW")

### lepidium
cat("lepidium\n")

CompareRobustMethods(data.to.comp = lepidium,
                             var.dose = "conc",
                             var.response = "weight")

### lettuce
cat("lettuce\n")

CompareRobustMethods(data.to.comp = lettuce,
                             var.dose = "conc",
                             var.response = "weight")

### mecter
cat("mecter\n")

CompareRobustMethods(data.to.comp = mecter,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr")

### M.bahia
cat("M.bahia\n")

CompareRobustMethods(data.to.comp = M.bahia,
                             var.dose = "conc",
                             var.response = "dryweight")

### nasturtium
cat("nasturtium\n")

CompareRobustMethods(data.to.comp = nasturtium,
                             var.dose = "conc",
                             var.response = "weight")

### O.mykiss
cat("O.mykiss\n")

CompareRobustMethods(data.to.comp = O.mykiss,
                             var.dose = "conc",
                             var.response = "weight")

### P.promelas
cat("P.promelas\n")

CompareRobustMethods(data.to.comp = P.promelas,
                             var.dose = "conc",
                             var.response = "dryweight")

### ryegrass
cat("ryegrass\n")

CompareRobustMethods(data.to.comp = ryegrass,
                             var.dose = "conc",
                             var.response = "rootl")

### S.alba
cat("S.alba\n")

CompareRobustMethods(data.to.comp = S.alba,
                             var.ref = "Herbicide",
                             var.dose = "Dose",
                             var.response = "DryMatter")

### S.capricornutum
cat("S.capricornutum\n")

CompareRobustMethods(data.to.comp = S.capricornutum,
                             var.dose = "conc",
                             var.response = "count")

### secalonic
cat("secalonic\n")

CompareRobustMethods(data.to.comp = secalonic,
                             var.dose = "dose",
                             var.response = "rootl")

### spinach
cat("spinach\n")

CompareRobustMethods(data.to.comp = spinach,
                             var.ref = "CURVE",
                             var.dose = "DOSE",
                             var.response = "SLOPE")

### terbuthylazin
cat("terbuthylazin\n")

CompareRobustMethods(data.to.comp = terbuthylazin,
                             var.dose = "dose",
                             var.response = "rgr")

### vinclozolin
cat("vinclozolin\n")

CompareRobustMethods(data.to.comp = vinclozolin,
                             var.ref = "exper",
                             var.dose = "conc",
                             var.response = "effect")

sink()

