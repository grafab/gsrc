context("Sample subsetting")

test_that("filt_samp returns correct raw_data object", {
  if(require(brassicaData)){
    data("raw_napus", package = "brassicaData")
    raw_napus <- filt_samp(raw_napus, check_raw(raw = raw_napus, plot = FALSE, thresh = 28000))
    expect_equal(class(raw_napus), "raw_data")
    expect_true(is.list(raw_napus))
    expect_equal(length(raw_napus$samples), 280)
  }
})