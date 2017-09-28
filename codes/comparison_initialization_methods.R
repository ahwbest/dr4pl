# -----------------------------------------------------------------------
### User-defined functions
#
CompareInitializationMethods <- function(data.input,
                                         var.dose,
                                         var.response,
                                         var.ref = NULL) {

  data.input <- na.omit(data.input)
  n.methods <- 2  # Number of initialization methods to be tried
  
  result.mat <- matrix(nrow = n.methods, ncol = 5)
  
  if(length(var.ref) == 0) {
    
    data.whole <- subset(x = data.input, select = c(var.dose, var.response))
    colnames(data.whole) <- c("Dose", "Response")
    
    data.part <- data.whole
    
    obj.drra.logistic <- drra(Response ~ Dose,
                              data = data.part,
                              method.init = "Logistic")
    result.mat[1, ] <- c(obj.drra.logistic$parameters, obj.drra.logistic$error.value)
    
    obj.drra.Mead <- drra(Response ~ Dose,
                          data = data.part,
                          method.init = "Mead")
    result.mat[2, ] <- c(obj.drra.Mead$parameters, obj.drra.Mead$error.value)
    
    row.names(result.mat) <- c("Logistic", "Mead")
    colnames(result.mat) <- c("Upper limit", "IC50", "Slope", "Lower limit", "Loss value")
    
    print(result.mat)
    cat("\n")
    
    # result.comparison <- result.comparison&
    #   (signif(result[1, 5], 3) == signif(result[2, 5], 3))
  } else {
    
    data.whole <- subset(x = data.input, select = c(var.ref, var.dose, var.response))
    colnames(data.whole) <- c("Ref", "Dose", "Response")
    data.whole$Ref <- as.factor(data.whole$Ref)
    data.whole <- subset(x = data.whole, subset = Ref != "999")
    data.whole <- droplevels(data.whole)
    
    levels.ref <- levels(data.whole$Ref)
    
    for(i in 1:length(levels.ref)) {
    
      data.part <- subset(x = data.whole, 
                          select = c(Dose, Response),
                          subset = Ref == levels.ref[i])

      obj.drra.logistic <- drra(Response ~ Dose,
                                data = data.part,
                                method.init = "Logistic")
      result.mat[1, ] <- c(obj.drra.logistic$parameters, obj.drra.logistic$error.value)
      
      obj.drra.Mead <- drra(Response ~ Dose,
                            data = data.part,
                            method.init = "Mead")
      result.mat[2, ] <- c(obj.drra.Mead$parameters, obj.drra.Mead$error.value)
      
      row.names(result.mat) <- c("Logistic", "Mead")
      colnames(result.mat) <- c("Upper limit", "IC50", "Slope", "Lower limit", "Loss value")
      
      print(result.mat)
      cat("\n")
    }
  }
}

sink(".\\tests\\testthat\\comparison_intialization_methods.Rout")

### acidiq
cat("acidiq\n")
CompareInitializationMethods(data.input = acidiq,
                       var.ref = "pct",
                       var.dose = "dose",
                       var.response = "rgr")

### algae
cat("algae\n")
CompareInitializationMethods(data.input = algae,
                       var.dose = "conc",
                       var.response = "vol")

### etmotc
cat("etmotc\n")
CompareInitializationMethods(data.input = etmotc,
                       var.ref = "pct1",
                       var.dose = "dose1",
                       var.response = "rgr1")

### G.aparine
cat("G.aparine\n")
CompareInitializationMethods(data.input = G.aparine,
                       var.dose = "dose",
                       var.response = "drymatter")

### glymet
cat("glylmet\n")
CompareInitializationMethods(data.input = glymet,
                       var.ref = "pct",
                       var.dose = "dose",
                       var.response = "rgr")

### heartrate
cat("heartrate\n")
CompareInitializationMethods(data.input = heartrate,
                       var.dose = "pressure",
                       var.response = "rate")

### leaflength
cat("leaflength\n")
CompareInitializationMethods(data.input = leaflength,
                       var.dose = "Dose",
                       var.response = "DW")

### lepidium
cat("lepidium\n")
CompareInitializationMethods(data.input = lepidium,
                       var.dose = "conc",
                       var.response = "weight")

### lettuce
cat("lettuce\n")
CompareInitializationMethods(data.input = lettuce,
                       var.dose = "conc",
                       var.response = "weight")

### mecter
cat("mecter\n")
CompareInitializationMethods(data.input = mecter,
                       var.ref = "pct",
                       var.dose = "dose",
                       var.response = "rgr")

### M.bahia
cat("M.bahia\n")
CompareInitializationMethods(data.input = M.bahia,
                       var.dose = "conc",
                       var.response = "dryweight")

### nasturtium
cat("nasturtium\n")
CompareInitializationMethods(data.input = nasturtium,
                       var.dose = "conc",
                       var.response = "weight")

### O.mykiss
cat("O.mykiss\n")
CompareInitializationMethods(data.input = O.mykiss,
                       var.dose = "conc",
                       var.response = "weight")

### P.promelas
cat("P.promelas\n")
CompareInitializationMethods(data.input = P.promelas,
                       var.dose = "conc",
                       var.response = "dryweight")

### ryegrass
cat("ryegrass\n")
CompareInitializationMethods(data.input = ryegrass,
                       var.dose = "conc",
                       var.response = "rootl")

### S.alba
cat("S.alba\n")
CompareInitializationMethods(data.input = S.alba,
                       var.ref = "Herbicide",
                       var.dose = "Dose",
                       var.response = "DryMatter")

### S.capricornutum
cat("S.capricornutum\n")
CompareInitializationMethods(data.input = S.capricornutum,
                       var.dose = "conc",
                       var.response = "count")

### secalonic
cat("secalonic\n")
CompareInitializationMethods(data.input = secalonic,
                       var.dose = "dose",
                       var.response = "rootl")

### spinach
cat("spinach\n")
CompareInitializationMethods(data.input = spinach,
                       var.ref = "CURVE",
                       var.dose = "DOSE",
                       var.response = "SLOPE")

### terbuthylazin
cat("terbuthylazin\n")
CompareInitializationMethods(data.input = terbuthylazin,
                       var.dose = "dose",
                       var.response = "rgr")

### vinclozolin
cat("vinclozolin\n")
CompareInitializationMethods(data.input = vinclozolin,
                       var.ref = "exper",
                       var.dose = "conc",
                       var.response = "effect")

sink()