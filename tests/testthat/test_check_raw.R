context("Check raw intensity values")

test_that("check_raw returns correct output", {
  if(require(brassicaData)){
    data("raw_napus", package = "brassicaData", envir = environment())
    to_filter <- check_raw(raw_napus, thresh = 28000, breaks = 20)
    expect_true(is.vector(to_filter))
    expect_equal(to_filter, c(1:18, 21, 22, 25, 26, 29, 30))
  }
})
