# vim:textwidth=80:expandtab:shiftwidth=2:softtabstop=2
library(whalestrike)
library(testthat)

context("whale properties")

test_that("parameters", {
  expect_equal(parameters()$lw, 13.7)
  expect_equal(parameters(lw = "from_species")$lw, 13.8)
  expect_equal(parameters()$l, c(0.025, 0.160, 1.120, 0.100))
  expect_equal(parameters(l = "from_species")$l, c(0.009, 0.163, 1.325, 0.143))
})

test_that("whaleMassFromLength() simple call as in app()", {
  expect_equal(whaleMassFromLength(13.7), 29993.89, tolerance=0.01)
  expect_equal(whaleMassFromLength(13.7, "Default"), 29993.89, tolerance=0.01)
  expect_equal(whaleMassFromLength(13.7, "N. Atl. Right Whale"), 29993.89, tolerance=0.01)
})

test_that("whaleMassFromLength and whaleLengthFromMass are inverses", {
  for (model in c("moore2005", "fortune2012atlantic", "fortune2012pacific")[1]) {
    L <- 1:20
    L <- 5
    M <- whaleMassFromLength(L, model = model)
    expect_equal(L, whaleLengthFromMass(M, model = model), tol = 1e-5)
  }
})

test_that("whale lethality index is self-consistent", {
  expect_equal(stressFromLethalityIndex(0.5), parameters()$logistic$tau50, tolerance = 0.01)
  expect_equal(lethalityIndexFromStress(parameters()$logistic$tau50), 0.5, tolerance = 0.01)
})

test_that("whale lethality index is unchanged (prevents ill-considered changes)", {
  expect_equal(stressFromLethalityIndex(0.5), 239883.3, scale = 1, tolerance = 0.1)
  expect_equal(lethalityIndexFromStress(239883.3), 0.5, scale = 1, tolerance = 0.000001)
})

test_that("strike produces same results as previously", {
  # Compare solution with results from package at the following state, which
  # seems to be the last commit before the final commit of the whale-collision
  # paper, so we can take it as a reference version.
  #  commit 84bde9dff8d3255e18fefbb5ff479484ce8239bb (HEAD)
  #  Author: dankelley <kelley.dan@gmail.com>
  #  Date:   Wed Jul 8 12:44:09 2020 -0300
  t <- seq(0, 0.8, length.out = 50)
  state <- list(xs = -2, vs = knot2mps(10), xw = 0, vw = 0)
  parms <- parameters()
  sol <- strike(t, state, parms)
  data(sol20200708)
  sol$parms$stressFromStrain <- NULL
  sol20200708$parms$stressFromStrain <- NULL
  expect_equal(sol, sol20200708)
})

test_that("shipMassFromLength", {
  expect_equal(shipMassFromLength("Tug", 50) / 1e3, 1920.648427)
})
