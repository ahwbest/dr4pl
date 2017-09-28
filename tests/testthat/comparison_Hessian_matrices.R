# ---------------------------------------------------------------------------------
### Load libraries
#
library(drc)

# ---------------------------------------------------------------------------------
### User-defined functions
#
CompareHessianMatrices <- function(data.input,
                                   var.dose,
                                   var.response,
                                   var.ref = NULL) {

  data.input <- na.omit(data.input)
  
  if(length(var.ref) == 0) {
    
    data.whole <- subset(x = data.input, select = c(var.dose, var.response))
    colnames(data.whole) <- c("Dose", "Response")
    
    data.part <- data.whole
    x <- data.part$Dose
    y <- data.part$Response
    
    ### Dose response model obtained by the logistic method
    obj.dr4pl.logistic <- dr4pl(Response ~ Dose,
                                data = data.part,
                                method.init = "logistic")
    hessian.cO <- obj.dr4pl.logistic$hessian
    hessian.dr4pl <- Hessian(obj.dr4pl.logistic$parameters, x, y)

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
      
      x <- data.part$Dose
      y <- data.part$Response
      
      ### Dose response model obtained by the logistic method
      obj.dr4pl.logistic <- dr4pl(Response ~ Dose,
                                  data = data.part,
                                  method.init = "logistic")
      hessian.cO <- obj.dr4pl.logistic$hessian
      hessian.dr4pl <- Hessian(obj.dr4pl.logistic$parameters, x, y)

    }
    
    n.comparisons <- length(levels.ref)
  }
  
  return(c(n.win.logistic, n.win.Mead, n.comparisons))
  
}

CompareHessianMatricesDRC <- function() {
  
  CompareHessianMatrices(data.input = acidiq,
                         var.ref = "pct",
                         var.dose = "dose",
                         var.response = "rgr")
  
  CompareHessianMatrices(data.input = algae,
                         var.dose = "conc",
                         var.response = "vol")
  
  CompareHessianMatrices(data.input = etmotc,
                         var.ref = "pct1",
                         var.dose = "dose1",
                         var.response = "rgr1")
  
  CompareHessianMatrices(data.input = G.aparine,
                         var.ref = "treatment",
                         var.dose = "dose",
                         var.response = "drymatter")
  
  CompareHessianMatrices(data.input = glymet,
                         var.ref = "pct",
                         var.dose = "dose",
                         var.response = "rgr")
  
  CompareHessianMatrices(data.input = heartrate,
                         var.dose = "pressure",
                         var.response = "rate")
  
  CompareHessianMatrices(data.input = leaflength,
                         var.dose = "Dose",
                         var.response = "DW")
  
  CompareHessianMatrices(data.input = lepidium,
                         var.dose = "conc",
                         var.response = "weight")
  
  CompareHessianMatrices(data.input = lettuce,
                         var.dose = "conc",
                         var.response = "weight")
  
  CompareHessianMatrices(data.input = mecter,
                         var.ref = "pct",
                         var.dose = "dose",
                         var.response = "rgr")
  
  CompareHessianMatrices(data.input = M.bahia,
                         var.dose = "conc",
                         var.response = "dryweight")
  
  CompareHessianMatrices(data.input = nasturtium,
                         var.dose = "conc",
                         var.response = "weight")
  
  CompareHessianMatrices(data.input = O.mykiss,
                         var.dose = "conc",
                         var.response = "weight")
  
  CompareHessianMatrices(data.input = P.promelas,
                         var.dose = "conc",
                         var.response = "dryweight")
  
  CompareHessianMatrices(data.input = ryegrass,
                         var.dose = "conc",
                         var.response = "rootl")
  
  CompareHessianMatrices(data.input = S.alba,
                         var.ref = "Herbicide",
                         var.dose = "Dose",
                         var.response = "DryMatter")
  
  CompareHessianMatrices(data.input = S.capricornutum,
                         var.dose = "conc",
                         var.response = "count")
  
  CompareHessianMatrices(data.input = secalonic,
                         var.dose = "dose",
                         var.response = "rootl")
  
  CompareHessianMatrices(data.input = spinach,
                         var.ref = "CURVE",
                         var.dose = "DOSE",
                         var.response = "SLOPE")
  
  CompareHessianMatrices(data.input = terbuthylazin,
                         var.dose = "dose",
                         var.response = "rgr")
  
  CompareHessianMatrices(data.input = vinclozolin,
                         var.ref = "exper",
                         var.dose = "conc",
                         var.response = "effect")
  
}

# ---------------------------------------------------------------------------------
### Test
#
context("Test whether the Hessian matrix computed by the function \`Hessian\` coincides
        with the Hessian matrix returned by the function \`constrOptim\` on the data
        sets of the R package \`drc\`.")

