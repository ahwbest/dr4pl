# -----------------------------------------------------------------------
### User-defined functions
#
CompareIntialization <- function(data.input,
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
    
    row.names(result.mat) <- c("drc", "drra")
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
      
      row.names(result.mat) <- c("drc", "drra")
      colnames(result.mat) <- c("Upper limit", "IC50", "Slope", "Lower limit", "Loss value")
      
      print(result.mat)
      cat("\n")
    }
  }
}

cat("acidiq\n")
CompareInitializationMethods(data.input = acidiq,
                             var.ref = "pct",
                             var.dose = "dose",
                             var.response = "rgr")

cat("algae\n")
CompareInitializationMethods(data.input = algae,
                             var.dose = "conc",
                             var.response = "vol")