
library(testthat)

context("Test whether all the error messages of the R package are correctly drawn
        when there are some errors in data or code.")

test_that("No error message is drawn.", {

  x <- 1:10
  theta <- c(100, -1, -1, 0)

  expect_error(MeanResponse(x, theta), "The IC50 parameter estimates become negative during the optimization process.")
})

# Compare the parameter estimates when the initial estimates are given by a user
# and computed by the Logistic method.
