# library(testthat)
# library(usethis)
# library(provoc)
#
# use_r("provoc_general")    # creates and opens R/foofy.R
# use_test("provoc_general") # creates and opens tests/testthat/test-blarg.R

test_that("importation of h5 files", {
  sp <- list("plip")
  expect_type(sp, type = "list")
})
