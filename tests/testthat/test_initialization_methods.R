# -----------------------------------------------------------------------
### User-defined functions
#
CompareInitializationMethods <- function(data.input,
                                         var.dose,
                                         var.response,
                                         var.ref = NULL) {

  data.input <- na.omit(data.input)
  n.methods <- 2  # Number of initialization methods to be tried
  n.comparisons <- 0
  n.win.logistic <- 0
  
  if(length(var.ref) == 0) {
    
    data.whole <- subset(x = data.input, select = c(var.dose, var.response))
    colnames(data.whole) <- c("Dose", "Response")
    
    data.part <- data.whole
    
    obj.drra.logistic <- drra(Response ~ Dose,
                              data = data.part,
                              method.init = "Logistic")
    loss.logistic <- obj.drra.logistic$error.value
    
    obj.drra.Mead <- drra(Response ~ Dose,
                          data = data.part,
                          method.init = "Mead")
    loss.mead <- obj.drra.Mead$error.value
    
    n.win.logistic <- n.win.logistic + as.numeric(loss.logistic<loss.mead)
    n.comparisons <- 1
    
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
      loss.logistic <- obj.drra.logistic$error.value
      
      obj.drra.Mead <- drra(Response ~ Dose,
                            data = data.part,
                            method.init = "Mead")
      loss.mead <- obj.drra.Mead$error.value
      
      n.win.logistic <- n.win.logistic + as.numeric(loss.logistic<loss.mead)
    }
    
    n.comparisons <- length(levels.ref)
  }
  
  return(c(n.win.logistic, n.comparisons))
}

TestInitializationMethods <- function() {
  
  results <- c(0, 0)
  
  results <- results + CompareInitializationMethods(data.input = acidiq,
                                                    var.ref = "pct",
                                                    var.dose = "dose",
                                                    var.response = "rgr")
  
  results <- results + CompareInitializationMethods(data.input = algae,
                                                    var.dose = "conc",
                                                    var.response = "vol")
  
  results <- results + CompareInitializationMethods(data.input = etmotc,
                                                    var.ref = "pct1",
                                                    var.dose = "dose1",
                                                    var.response = "rgr1")
  
  results <- results + CompareInitializationMethods(data.input = G.aparine,
                                                    var.dose = "dose",
                                                    var.response = "drymatter")
  
  results <- results + CompareInitializationMethods(data.input = glymet,
                                                    var.ref = "pct",
                                                    var.dose = "dose",
                                                    var.response = "rgr")
  
  results <- results + CompareInitializationMethods(data.input = heartrate,
                                                    var.dose = "pressure",
                                                    var.response = "rate")
  
  results <- results + CompareInitializationMethods(data.input = leaflength,
                                                    var.dose = "Dose",
                                                    var.response = "DW")
  
  results <- results + CompareInitializationMethods(data.input = lepidium,
                                                    var.dose = "conc",
                                                    var.response = "weight")
  
  results <- results + CompareInitializationMethods(data.input = lettuce,
                                                    var.dose = "conc",
                                                    var.response = "weight")
  
  results <- results + CompareInitializationMethods(data.input = mecter,
                                                    var.ref = "pct",
                                                    var.dose = "dose",
                                                    var.response = "rgr")
  
  results <- results + CompareInitializationMethods(data.input = M.bahia,
                                                    var.dose = "conc",
                                                    var.response = "dryweight")
  
  results <- results + CompareInitializationMethods(data.input = nasturtium,
                                                    var.dose = "conc",
                                                    var.response = "weight")
  
  results <- results + CompareInitializationMethods(data.input = O.mykiss,
                                                    var.dose = "conc",
                                                    var.response = "weight")
  
  results <- results + CompareInitializationMethods(data.input = P.promelas,
                                                    var.dose = "conc",
                                                    var.response = "dryweight")
  
  results <- results + CompareInitializationMethods(data.input = ryegrass,
                                                    var.dose = "conc",
                                                    var.response = "rootl")
  
  results <- results + CompareInitializationMethods(data.input = S.alba,
                                                    var.ref = "Herbicide",
                                                    var.dose = "Dose",
                                                    var.response = "DryMatter")
  
  results <- results + CompareInitializationMethods(data.input = S.capricornutum,
                                                    var.dose = "conc",
                                                    var.response = "count")
  
  results <- results + CompareInitializationMethods(data.input = secalonic,
                                                    var.dose = "dose",
                                                    var.response = "rootl")
  
  results <- results + CompareInitializationMethods(data.input = spinach,
                                                    var.ref = "CURVE",
                                                    var.dose = "DOSE",
                                                    var.response = "SLOPE")
  
  results <- results + CompareInitializationMethods(data.input = terbuthylazin,
                                                    var.dose = "dose",
                                                    var.response = "rgr")
  
  results <- results + CompareInitializationMethods(data.input = vinclozolin,
                                                    var.ref = "exper",
                                                    var.dose = "conc",
                                                    var.response = "effect")
  
  return(results[1])
}

test_that("The number of times that the logistic method outperforms Mead's method
           on the `drc' data sets should remain the same.", {
  
  expect_equal(TestInitializationMethods(), 3)
           })