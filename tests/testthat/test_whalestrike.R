## vim:textwidth=80:expandtab:shiftwidth=2:softtabstop=2
library(whalestrike)
library(testthat)

context("whale properties")

test_that("whaleMassFromLength and whaleLengthFromMass are inverses", {
          for (model in c("moore2005", "fortune2012atlantic", "fortune2012pacific")[1]) {
            L <- 1:20
            L <- 5
            M <- whaleMassFromLength(L, model=model)
            expect_equal(L, whaleLengthFromMass(M, model=model), tol=1e-5)
          }
})

