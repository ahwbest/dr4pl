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

test_that("Error messages are correctly drawn.", {

  x <- 1:10  # Doses
  theta <- c(100, -1, -1, 0)  # Parameters of the 4PL model

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
  
  ### The initialization and robust estimation methods should be correctly specified.
  expect_error(FindInitialParms(x, y, "logistic", "abc"))
  expect_error(FindInitialParms(x, y, "abc", "absolute"))
  
  ### The title text of the plot function should be of the character type.
  expect_error(plot(dr4pl.test, text.title = 123.45))
})

test_that("Test whether the conventional R packages for dose-response modelling such
          as `drc' and `nplr' fail in error data sets.", {
            
            expect_error(drm(Response ~ Dose, data = drc_error_1,
                             fct = LL.4()))
            expect_error(drm(Response ~ Dose, data = drc_error_2,
                             fct = LL.4()))
            expect_error(drm(Response ~ Dose, data = drc_error_3,
                             fct = LL.4()))
            expect_error(drm(Response ~ Dose, data = drc_error_4,
                             fct = LL.4()))
            
          })