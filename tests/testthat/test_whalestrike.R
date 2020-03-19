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

test_that("whale lethality index is self-consistent", {
          expect_equal(stressFromLethalityIndex(0.5), parameters()$logistic$tau50, tolerance=0.01)
          expect_equal(lethalityIndexFromStress(parameters()$logistic$tau50), 0.5, tolerance=0.01)
})

test_that("whale lethality index is unchanged (prevents ill-considered changes)", {
          expect_equal(stressFromLethalityIndex(0.5), 239883.3, scale=1, tolerance=0.1)
          expect_equal(lethalityIndexFromStress(239883.3), 0.5, scale=1, tolerance=0.000001)
})

