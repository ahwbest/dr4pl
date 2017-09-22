# ----------------------------------------------------------------------------
### Load data
#
n.error.cases <- 4  # Number of error cases.

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

data.error.list <- vector("list", n.error.cases)

for(i in 1:n.error.cases) {
  
  file.name <- paste("drc_error_", i, ".csv", sep = "")
  data.error.list[[i]] <- read.csv(file.name)
}

# ---------------------------------------------------------------------------
### Test
#
for(i in 1:n.error.cases) {
  
  data.error <- data.error.list[[i]]
  
  obj.dr4pl.logistic <- dr4pl(Response ~ Dose,
                              data = data.error,
                              method.init = "logistic")
  parms.logistic <- obj.dr4pl.logistic$parameters
  loss.logistic <- signif(obj.dr4pl.logistic$error.value, 4)
    
  obj.dr4pl.Mead <- dr4pl(Response ~ Dose,
                          data = data.error,
                          method.init = "Mead")
  parms.Mead <- obj.dr4pl.Mead$parameters
  loss.Mead <- signif(obj.dr4pl.Mead$error.value, 4)
  
  cat(paste(c(loss.logistic, loss.Mead), "\n"))
}



### Dose response model obtained by the Mead method


