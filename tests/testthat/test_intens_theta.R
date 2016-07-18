context("Intensity and Theta calculation")

test_that("intens_theta returns correct output", {
  if(require(brassicaData)){
    data("raw_napus", package = "brassicaData")
    raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:10)])
    raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:10)])
    dat <- intens_theta(raw_napus)
    expect_equal(class(dat), "norm_data")
    expect_true(is.list(dat))
    expect_true(is.matrix(dat$theta))
    expect_true(is.matrix(dat$intensity))
    expect_equal(dim(dat$theta), dim(dat$intensity))
  }
})