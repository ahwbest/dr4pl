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

  ### Check the errors of `dr4pl' are correctly drawn
  theta.init.trial <- c(100, -1, -1, 0)
  
  data.test <- data.frame(x = c(0.0001, 0.001, 0.01, 0.1, 1),
                          y = c(10, 9, 5, 1, 0))
  
  expect_error(dr4pl(y ~ x,
                     data = data.test,
                     init.parm = theta.init.trial,
                     method.init = "logistic"),
               "The IC50 parameter should be positive.")
  expect_error(dr4pl(y ~ x,
                     data = data.test,
                     method.init = "abc"),
               "The initialization method name should be one of \"logistic\" and \"Mead\".")
  expect_error(dr4pl(y ~ x,
                     data = data.test,
                     trend = "abc"),
               "The type of the \"trend\" parameter should be one of \"auto\", \"decreasing\" and \"increasing\".")
  data.test$x[1] <- -1
  expect_error(dr4pl(y ~ x,
                     data = data.test),
               "Dose levels should be nonnegative.")
  data.test$x[1] <- 1
  
  ### Check the errors of `FindInitialParms' are correctly drawn
  expect_error(FindInitialParms(x, y, "logistic", "abc"))
  expect_error(FindInitialParms(x, y, "abc", "absolute"))
  
  ### Check the errors of `plot.dr4pl` are correctly drawn
  expect_error(plot(dr4pl(y ~ x,
                          data = data.test,
                          method.init = "logistic"), text.title = 143.45),
               "Title text should be characters.")
})