# --------------------------------------------------------------------------------
### Load libraries
#
library(testthat)

# -------------------------------------------------------------------------------
### Test
#
context("Test whether all the error messages of the R package are correctly drawn
        when there are some errors in data or code.")

test_that("No error message is drawn.", {

  x <- 1:10
  theta <- c(100, -1, -1, 0)

  expect_error(MeanResponse(x, theta), "The IC50 parameter estimates become negative during the optimization process.")
  
  theta.init.trial <- c(100, -1, -1, 0)
  
  data.test <- data.frame(x = c(0.0001, 0.001, 0.01, 0.1, 1),
                          y = c(10, 9, 5, 1, 0))
  
  expect_error(dr4pl(y ~ x,
                     data = data.test,
                     init.parm = theta.init.trial,
                     method.init = "logistic"),
               "Initial parameter values are not in the interior of the feasible region.")
  
  dr4pl.test <- dr4pl(y ~ x,
                      data = data.test,
                      method.init = "logistic")
  
  expect_error(plot(dr4pl.test, text.title = 123.45),
               "Initial parameter values are not in the interior of the feasible region.")
})