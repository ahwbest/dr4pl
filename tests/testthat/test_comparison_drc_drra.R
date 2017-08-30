# -----------------------------------------------------------------------------
### Load libraries and source code
#
library(drc)
library(ggplot2)
library(RColorBrewer)
library(drra)
library(testthat)

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
  #data(list=data.to.comp, package="drc")
  data.to.comp <- na.omit(data.to.comp)

  result.comparison <- TRUE

  # There is only one drug
  if(length(var.ref) == 0) {

    data.new <- subset(x = data.to.comp, select = c(var.dose, var.response))
    colnames(data.new) <- c("Dose", "Response")
    n <- nrow(data.new)  # Number of observations

    # drc
    obj.drc <- drc::drm(Response ~ Dose,
                     data = data.new,
                     fct = drc::LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
                                method = "1"),
                     control = drc::drmc(method = "Nelder-Mead"))

    result <- coef(obj.drc)
    result <- result[c(3, 4, 1, 2)]
    result[3] <- result[3]

    result <- c(result, obj.drc$fit$value/n)

    # drra
    obj.drra <- drra(Response ~ Dose,
                     data = data.new)

    parm.drra <- coef(obj.drra)

    result <- rbind(result, c(parm.drra, obj.drra$error.value))

    row.names(result) <- c("drc", "drra")
    colnames(result) <- c("Lower limit", "IC50", "Slope", "Upper limit", "Loss value")

    print(result)
    cat("\n")

    result.comparison <- result.comparison&
                         (signif(result[1, 5], 3) == signif(result[2, 5], 3))

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

  } else {  # There are multiple drugs

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
      drc.ctrl <- drc::drmc(method = "Nelder-Mead")

      obj.drc <- drc::drm(Response ~ Dose,
                       data = data.each,
                       fct = drc::LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50")),
                       control = drc.ctrl)

      result <- coef(obj.drc)
      result <- result[ c(3, 4, 1, 2)]
      result[3] <- -result[3]

      result <- c(result, obj.drc$fit$value/n)

      # drra
      obj.drra <- drra(Response ~ Dose,
                       data = data.each)

      parm.drra <- coef(obj.drra)

      result <- rbind(result, c(parm.drra, obj.drra$error.value))

      row.names(result) <- c("drc", "drra")
      colnames(result) <- c("Left limit", "IC50", "Slope", "Right limit", "Loss value")

      print(result)
      cat("\n")

      result.comparison <- result.comparison&
                           (signif(result[1, 5], 3) == signif(result[2, 5], 3))
    }
  }

  return(result.comparison)
}

# -------------------------------------------------------------------------------
### Compare the parameter estimates
#
### acidiq
test_that("The output of drc and drra should be the same.", {

cat("acidiq\n")
expect_true(CompareInitializationMethods(data.to.comp = acidiq,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr"))

### algae
cat("algae\n")
# Plot fitted curves
expect_true(CompareInitializationMethods(data.to.comp = algae,
                             var.dose = "conc",
                             var.response = "vol"))

### etmotc
cat("etmotc\n")
expect_true(CompareInitializationMethods(data.to.comp = etmotc,
                             var.ref = "pct1",
                             var.dose = "dose1",
                             var.response = "rgr1"))

### G.aparine
cat("G.aparine\n")
expect_true(CompareInitializationMethods(data.to.comp = G.aparine,
                             var.ref = "treatment",
                             var.dose = "dose",
                             var.response = "drymatter"))

### glymet
cat("glylmet\n")
expect_true(CompareInitializationMethods(data.to.comp = glymet,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr"))

### heartrate
cat("heartrate\n")
expect_true(CompareInitializationMethods(data.to.comp = heartrate,
                             var.dose = "pressure",
                             var.response = "rate"))

### leaflength
cat("leaflength\n")
expect_true(CompareInitializationMethods(data.to.comp = leaflength,
                             var.dose = "Dose",
                             var.response = "DW"))

### lepidium
cat("lepidium\n")
expect_true(CompareInitializationMethods(data.to.comp = lepidium,
                             var.dose = "conc",
                             var.response = "weight"))

### lettuce
cat("lettuce\n")
expect_true(CompareInitializationMethods(data.to.comp = lettuce,
                             var.dose = "conc",
                             var.response = "weight"))

### mecter
cat("mecter\n")
expect_true(CompareInitializationMethods(data.to.comp = mecter,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr"))

### M.bahia
cat("M.bahia\n")
expect_true(CompareInitializationMethods(data.to.comp = M.bahia,
                             var.dose = "conc",
                             var.response = "dryweight"))

### nasturtium
cat("nasturtium\n")
expect_true(CompareInitializationMethods(data.to.comp = nasturtium,
                             var.dose = "conc",
                             var.response = "weight"))

### O.mykiss
cat("O.mykiss\n")
expect_true(CompareInitializationMethods(data.to.comp = O.mykiss,
                             var.dose = "conc",
                             var.response = "weight"))

### P.promelas
cat("P.promelas\n")
expect_true(CompareInitializationMethods(data.to.comp = P.promelas,
                             var.dose = "conc",
                             var.response = "dryweight"))

### ryegrass
cat("ryegrass\n")
expect_true(CompareInitializationMethods(data.to.comp = ryegrass,
                             var.dose = "conc",
                             var.response = "rootl"))

### S.alba
cat("S.alba\n")
expect_true(CompareInitializationMethods(data.to.comp = S.alba,
                             var.ref = "Herbicide",
                             var.dose = "Dose",
                             var.response = "DryMatter"))

### S.capricornutum
cat("S.capricornutum\n")
expect_true(CompareInitializationMethods(data.to.comp = S.capricornutum,
                             var.dose = "conc",
                             var.response = "count"))

### secalonic
cat("secalonic\n")
expect_true(CompareInitializationMethods(data.to.comp = secalonic,
                             var.dose = "dose",
                             var.response = "rootl"))

### spinach
cat("spinach\n")
expect_true(CompareInitializationMethods(data.to.comp = spinach,
                             var.ref = "CURVE",
                             var.dose = "DOSE",
                             var.response = "SLOPE"))

### terbuthylazin
cat("terbuthylazin\n")
expect_true(CompareInitializationMethods(data.to.comp = terbuthylazin,
                             var.dose = "dose",
                             var.response = "rgr"))

### vinclozolin
cat("vinclozolin\n")
expect_true(CompareInitializationMethods(data.to.comp = vinclozolin,
                             var.ref = "exper",
                             var.dose = "conc",
                             var.response = "effect"))
})
