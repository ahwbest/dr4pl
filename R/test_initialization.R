# -----------------------------------------------------------------------------
### Load libraries and source code
#
library(drc)
source(".\\R\\initialization.R")
source(".\\R\\main.R")

CompareInitializationMethods <- function(data.to.comp,
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

  result <- c()

  if(length(var.ref) == 0) {

    data.new <- subset(x = data.to.comp, select = c(var.dose, var.response))
    colnames(data.new) <- c("Dose", "Response")
    n <- nrow(data.new)  # Number of observations

    # drc
    obj.drc.1 <- drm(Response ~ Dose,
                     data = data.new,
                     fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
                                method = "1"),
                     control = drmc(method = "Nelder-Mead"))

    obj.drc.2 <- drm(Response ~ Dose,
                     data = data.new,
                     fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
                                method = "2"),
                     control = drmc(method = "Nelder-Mead"))

    obj.drc.3 <- drm(Response ~ Dose,
                     data = data.new,
                     fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
                                method = "3"),
                     control = drmc(method = "Nelder-Mead"))

    result <- rbind(coef(obj.drc.1), coef(obj.drc.2), coef(obj.drc.3))
    result <- result[, c(3, 4, 1, 2)]
    result[, 3] <- -result[, 3]
    result <- cbind(result, c(obj.drc.1$fit$value, obj.drc.2$fit$value, obj.drc.3$fit$value)/n)

    # drra
    obj.drra <- drra(Response ~ Dose,
                     data = data.new)

    parm.drra <- coef(obj.drra)

    result <- rbind(result, c(parm.drra, obj.drra$error.value))

    row.names(result) <- c("drc 1", "drc 2", "drc 3", "drra")
    colnames(result) <- c("Lower limit", "IC50", "Slope", "Upper limit", "Loss value")

    print(result)
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
    #   ComparisonPlot(dose = data.new$Dose, response = data.new$Response, parm.mat)
    # }

  } else {

    data.new <- subset(x = data.to.comp, select = c(var.ref, var.dose, var.response))
    colnames(data.new) <- c("Ref", "Dose", "Response")
    data.new$Ref <- as.factor(data.new$Ref)
    data.new <- subset(x = data.new, subset = Ref != "999")
    data.new <- droplevels(data.new)

    level.ref <- levels(data.new$Ref)

    for(i in 1:length(level.ref)) {
      data.each <- subset(x = data.new, subset = Ref == level.ref[i])
      n <- nrow(data.each)  # Number of observations

      cat(var.ref, "=", level.ref[i], "\n")

      # drc
      obj.drc.1 <- drm(Response ~ Dose,
                       data = data.each,
                       fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
                                  method = "1"),
                       control = drmc(method = "Nelder-Mead"))

      obj.drc.2 <- drm(Response ~ Dose,
                       data = data.each,
                       fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
                                  method = "2"),
                       control = drmc(method = "Nelder-Mead"))

      obj.drc.3 <- drm(Response ~ Dose,
                       data = data.each,
                       fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
                                  method = "3"),
                       control = drmc(method = "Nelder-Mead"))

      result <- rbind(coef(obj.drc.1), coef(obj.drc.2), coef(obj.drc.3))
      result <- result[, c(3, 4, 1, 2)]
      result[, 3] <- -result[, 3]
      result <- cbind(result, c(obj.drc.1$fit$value, obj.drc.2$fit$value, obj.drc.3$fit$value)/n)

      # drra
      obj.drra <- drra(Response ~ Dose,
                       data = data.each)

      parm.drra <- coef(obj.drra)

      result <- rbind(result, c(parm.drra, obj.drra$error.value))

      row.names(result) <- c("drc 1", "drc 2", "drc 3", "drra")
      colnames(result) <- c("Left limit", "IC50", "Slope", "Right limit", "Loss value")

      print(result)
      cat("\n")

    }
  }
}

# -------------------------------------------------------------------------------
### Compare the parameter estimates
#
sink(file = "parameter_comparison.txt")

### acidiq
cat("acidiq\n")

CompareInitializationMethods(data.to.comp = acidiq,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr")

### algae
cat("algae\n")

# Plot fitted curves
CompareInitializationMethods(data.to.comp = algae,
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
                             var.ref = "treatment",
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
