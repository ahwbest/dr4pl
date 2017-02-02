# -----------------------------------------------------------------------------
### Load libraries and source code
#
library(drc)
library(ggplot2)
library(RColorBrewer)
library(drra)

ComparisonPlot <- function(dose, response, parm.matr) {
# Draw curves of different 4PL models
#
# Args:
#   dose: Vector of dose values
#   response: Vector of response values
#   parm.matr: Matrix of parameter values. Each row corresponds to paramters of
#              of each model. There are 4 columns corresponding to 4 parameters.

  if(length(dose) != length(response)) {
    stop("The numbers of does and response values should be the same.")
  }

  n.curves <- nrow(parm.matr)  # Number of curves to be plotted
  data.plot <- data.frame(Dose = dose, Response = response)
  color.vec <- brewer.pal(n.curves, "Dark2")  # Obtain colors from brewer.pal

  a <- ggplot(aes(x = dose, y = response), data = data.plot)

  for(i in 1:n.curves) {

    a <- a + stat_function(fun = MeanResponse,
                           args = list(theta = parm.matr[i, ]),
                           colour = color.vec[i],
                           size = 1.2)
  }

  a <- a + geom_point(size = I(5), alpha = I(0.8))

  a <- a + labs(title = "Dose response curve")
  a <- a + xlab("Dose") + ylab("Response")

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

  plot(a)
}

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

  result <- c()

  # There is only one drug
  if(length(var.ref) == 0) {

    data.new <- subset(x = data.to.comp, select = c(var.dose, var.response))
    colnames(data.new) <- c("Dose", "Response")
    n <- nrow(data.new)  # Number of observations

    # drc
    obj.drc.1 <- drc::drm(Response ~ Dose,
                     data = data.new,
                     fct = drc::LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
                                method = "1"),
                     control = drc::drmc(method = "Nelder-Mead"))

    result <- rbind(coef(obj.drc.1))
    #result <- rbind(coef(obj.drc.1), coef(obj.drc.2), coef(obj.drc.3))
    #result <- result[, c(3, 4, 1, 2)]
    result <- -result
    result <- cbind(result, c(obj.drc.1$fit$value)/n)

    # drra
    obj.drra <- drra(Response ~ Dose,
                     data = data.new)

    parm.drra <- coef(obj.drra)

    result <- rbind(result, c(parm.drra, obj.drra$error.value))

    row.names(result) <- c("drc 1", "drra")
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
      drc.ctrl <- drc::drmc(method = "Nelder-Mead", trace = TRUE)

      sink("output_drc.txt")

      obj.drc.1 <- drc::drm(Response ~ Dose,
                       data = data.each,
                       fct = drc::LL.4(names = c("Slope", "Lower limit", "Upper limit", "IC50"),
                                  method = "1"),
                       control = drc.ctrl)

      sink()

      result <- rbind(coef(obj.drc.1))
      #result <- result[, c(3, 4, 1, 2)]
      result <- -result
      result <- cbind(result, c(obj.drc.1$fit$value)/n)

      # drra
      sink("output_drra.txt")

      obj.drra <- drra(Response ~ Dose,
                       data = data.each)

      sink()

      parm.drra <- coef(obj.drra)

      result <- rbind(result, c(parm.drra, obj.drra$error.value))

      row.names(result) <- c("drc 1", "drra")
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

files_list<-c("G.aparine", "H.virescens", "M.bahia", "O.mykiss", "P.promelas", "RScompetition",
              "S.alba", "S.capricornutum", "acidiq", "algae", "auxins", "chickweed", "chickweed0",
              "daphnids", "decontaminants", "deguelin", "earthworms", "etmotc", "finney71",
              "germination", "glymet", "heartrate", "leaflength", "lepidium", "lettuce", "mecter",
              "metals", "methionine", "nasturtium", "ryegrass", "secalonic", "selenium", "spinach",
              "terbuthylazin", "vinclozolin"
)
lapply(files_list, function(filename){
  data(list=filename, package = "drc")
})

### acidiq
cat("acidiq\n")

CompareInitializationMethods(data.to.comp = acidiq,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr")
#rm(acidiq)

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
