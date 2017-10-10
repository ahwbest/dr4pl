## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(dr4pl)
library(ggplot2)
library(drc)

## ------------------------------------------------------------------------
a <- tryCatch({
  drc::drm(Response~Dose, data = drc_error_1, fct = LL.4())
}, 
warning = function(war) {
  # warning handler picks up where error was generated
  print(paste(war))
},
error = function(err) {
  # error handler picks up where error was generated
  print(paste(err))
})

## ------------------------------------------------------------------------
a <-tryCatch({
  dr4pl(Response~Dose, data = drc_error_1, method.robust = "Tukey")
},
warning = function(war) {
    # warning handler picks up where error was generated
    print(paste(war))
},
error = function(err) {
  # error handler picks up where error was generated
  print(paste(err))
})
plot(a, text.title = "Error plot #1", indices.outlier = c(102))

## ------------------------------------------------------------------------
a <- tryCatch({
  drc::drm(Response~Dose, data = drc_error_2, fct = LL.4())
}, 
warning = function(war) {
  # warning handler picks up where error was generated
  print(paste(war))
},
error = function(err) {
  # error handler picks up where error was generated
  print(paste(err))
})

## ------------------------------------------------------------------------
a <-tryCatch({
  dr4pl(Response~Dose, data = drc_error_2, method.init = "Mead", method.robust = "Huber" )
},
warning = function(war) {
    # warning handler picks up where error was generated
    print(paste(war))
},
error = function(err) {
  # error handler picks up where error was generated
  print(paste(err))
})
b <- plot(a, breaks.x = c(0.00135, 0.0135, 0.135, 1.35, 13.5), text.title = "Error plot #2", indices.outlier = c(2,8) )
b

## ------------------------------------------------------------------------
a <- tryCatch({
  drc::drm(Response~Dose, data = drc_error_3, fct = LL.4())
}, 
warning = function(war) {
  # warning handler picks up where error was generated
  print(paste(war))
},
error = function(err) {
  # error handler picks up where error was generated
  print(paste(err))
})


## ------------------------------------------------------------------------
a <-tryCatch({
  dr4pl(Response~Dose, data = drc_error_3, method.init = "Mead", method.robust = "absolute" )
},
warning = function(war) {
    # warning handler picks up where error was generated
    print(paste(war))
},
error = function(err) {
  # error handler picks up where error was generated
  print(paste(err))
})
plot(a, indices.outlier = c(90, 101), text.title = "Error plot #3")

## ------------------------------------------------------------------------
a <- tryCatch({
  drc::drm(Response~Dose, data = drc_error_4, fct = LL.4())
}, 
warning = function(war) {
  # warning handler picks up where error was generated
  print(paste(war))
},
error = function(err) {
  # error handler picks up where error was generated
  print(paste(err))
})

## ------------------------------------------------------------------------
a <-tryCatch({
  dr4pl(Response~Dose, data = drc_error_4, method.init = "Mead", method.robust = "absolute" )
},
warning = function(war) {
    # warning handler picks up where error was generated
    print(paste(war))
},
error = function(err) {
  # error handler picks up where error was generated
  print(paste(err))
})
plot(a, text.title = "Error plot #4", indices.outlier = c(1,100))

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_1,  method.init = "Mead")
plot(a, text.title = "Sample plot #1")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_2,  method.init = "Mead")
plot(a, text.title = "Sample plot #2")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_3, method.init = "Mead")
plot(a, text.title = "Sample plot #3")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_4, method.init = "Mead")
plot(a, text.title = "Sample plot #4")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_5, method.init = "Mead")
plot(a, text.title = "Sample plot #5")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_6, method.init = "Mead")
plot(a, text.title = "Sample plot #6")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_7, method.init = "Mead")
plot(a, text.title = "Sample plot #7")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_8, method.init = "Mead")
plot(a, text.title = "Sample plot #8")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_9, method.init = "Mead")
plot(a, text.title = "Sample plot #9")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_10, method.init = "Mead")
plot(a, text.title = "Sample plot #10")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_11, method.init = "Mead")
plot(a, text.title = "Sample plot #11")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_12, method.init = "Mead")
plot(a, text.title = "Sample plot #12")

## ------------------------------------------------------------------------
a <- dr4pl(Response~Dose, data = sample_data_13)
plot(a, text.title = "Sample plot #13")

