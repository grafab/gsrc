context("SNP subsetting")

test_that("filt_snp returns correct data type", {
  if(require(brassicaData)){
    data("raw_napus", package = "brassicaData")
    raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:10)])
    expect_equal(class(raw_napus), "raw_data")
    expect_equal(mode(raw_napus), "list")
    expect_equal(length(raw_napus$samples), 10)
    dat <- intens_theta(raw_napus)
    dat <- remove_suffix(dat, "_Grn")
    dat <- geno_baf_rratio(dat, delthresh = 11)
    dat <- filt_snps(dat, dat$snps[is.na(rowMeans(dat$baf, na.rm = TRUE))])
    expect_equal(class(dat), "norm_data")
    expect_true(is.list(dat))
  }
})