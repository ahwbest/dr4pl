# -----------------------------------------------------------------------------
### Load libraries
#
library(drc)
library(testthat)

# -----------------------------------------------------------------------------
### User-defined functions
#
Fit4PLVariousOptions <- function(data.input,
                                 var.dose,
                                 var.response,
                                 var.ref = NULL) {

  types.method.init <- c("logistic", "Mead")
  types.method.optim <- c("Nelder-Mead", "BFGS", "CG", "SANN")
  
  data.input <- na.omit(data.input)  # Omit NA values from data

  if(length(var.ref) == 0) {
    
    data.whole <- subset(x = data.input, select = c(var.dose, var.response))
    colnames(data.whole) <- c("Dose", "Response")
    
    data.part <- data.whole
    
    for(method.init in types.method.init) {
      
      for(method.optim in types.method.optim) {
        
        dr4pl(Response ~ Dose,
              data = data.part,
              method.init = method.init,
              method.optim = method.optim)
      }
    }

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

      for(method.init in types.method.init) {
        
        for(method.optim in types.method.optim) {
          
          dr4pl(Response ~ Dose,
                data = data.part,
                method.init = method.init,
                method.optim = method.optim)
        }
      }
    }
  }
}

TestNormalCases <- function() {
  
  Fit4PLVariousOptions(data.input = acidiq,
                       var.ref = "pct",
                       var.dose = "dose",
                       var.response = "rgr")
  
  Fit4PLVariousOptions(data.input = algae,
                       var.dose = "conc",
                       var.response = "vol")
  
  Fit4PLVariousOptions(data.input = etmotc,
                       var.ref = "pct1",
                       var.dose = "dose1",
                       var.response = "rgr1")
  
  Fit4PLVariousOptions(data.input = G.aparine,
                       var.dose = "dose",
                       var.response = "drymatter")
  
  Fit4PLVariousOptions(data.input = glymet,
                       var.ref = "pct",
                       var.dose = "dose",
                       var.response = "rgr")
  
  Fit4PLVariousOptions(data.input = heartrate,
                       var.dose = "pressure",
                       var.response = "rate")
  
  Fit4PLVariousOptions(data.input = leaflength,
                       var.dose = "Dose",
                       var.response = "DW")
  
  Fit4PLVariousOptions(data.input = lepidium,
                       var.dose = "conc",
                       var.response = "weight")
  
  Fit4PLVariousOptions(data.input = lettuce,
                       var.dose = "conc",
                       var.response = "weight")
  
  Fit4PLVariousOptions(data.input = mecter,
                       var.ref = "pct",
                       var.dose = "dose",
                       var.response = "rgr")
  
  Fit4PLVariousOptions(data.input = M.bahia,
                       var.dose = "conc",
                       var.response = "dryweight")
  
  Fit4PLVariousOptions(data.input = nasturtium,
                       var.dose = "conc",
                       var.response = "weight")
  
  Fit4PLVariousOptions(data.input = O.mykiss,
                       var.dose = "conc",
                       var.response = "weight")
  
  Fit4PLVariousOptions(data.input = P.promelas,
                       var.dose = "conc",
                       var.response = "dryweight")
  
  Fit4PLVariousOptions(data.input = ryegrass,
                       var.dose = "conc",
                       var.response = "rootl")
  
  Fit4PLVariousOptions(data.input = S.alba,
                       var.ref = "Herbicide",
                       var.dose = "Dose",
                       var.response = "DryMatter")
  
  Fit4PLVariousOptions(data.input = S.capricornutum,
                       var.dose = "conc",
                       var.response = "count")
  
  Fit4PLVariousOptions(data.input = secalonic,
                       var.dose = "dose",
                       var.response = "rootl")
  
  Fit4PLVariousOptions(data.input = spinach,
                       var.ref = "CURVE",
                       var.dose = "DOSE",
                       var.response = "SLOPE")
  
  Fit4PLVariousOptions(data.input = terbuthylazin,
                       var.dose = "dose",
                       var.response = "rgr")
  
  Fit4PLVariousOptions(data.input = vinclozolin,
                       var.ref = "exper",
                       var.dose = "conc",
                       var.response = "effect")
}

# --------------------------------------------------------------------------------
### Test
#
context("Test whether running `dr4pl' on the data sets of `drc' using different
        parameter settings does not draw any error")

test_that("Running `dr4pl' on the data sets of `drc' using different
           parameter settings should not draw any error.", {
  
  expect_error(TestNormalCases(), NA)
})