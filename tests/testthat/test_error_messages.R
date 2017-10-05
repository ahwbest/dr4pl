# --------------------------------------------------------------------------------
### Load libraries
#
library(drc)
library(testthat)

# -------------------------------------------------------------------------------
### Test
#
context("Test whether all the error messages of the R package are correctly drawn
        when there are some errors in data or code")

test_that("Error messages are correctly drawn", {

  x <- 1:10  # Doses
  theta <- c(100, -1, -1, 0)  # Parameters of the 4PL model

  theta.init.trial <- c(100, -1, -1, 0)
  
  data.test <- data.frame(x = c(0.0001, 0.001, 0.01, 0.1, 1),
                          y = c(10, 9, 5, 1, 0))
  
  expect_error(dr4pl(y ~ x,
                     data = data.test,
                     init.parm = theta.init.trial,
                     method.init = "logistic"),
               "The IC50 parameter should be positive.")
  
  dr4pl.test <- dr4pl(y ~ x,
                      data = data.test,
                      method.init = "logistic")
  
  ### The initialization and robust estimation methods should be correctly specified.
  expect_error(FindInitialParms(x, y, "logistic", "abc"))
  expect_error(FindInitialParms(x, y, "abc", "absolute"))
  
  ### The title text of the plot function should be of the character type.
  expect_error(plot(dr4pl.test, text.title = 143.45))
})