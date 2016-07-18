context("Genotyping, B-allele frequency and Log R ratio calculation")

test_that("geno_baf_rratio returns correct output", {
  if(require(brassicaData)){
    data("raw_napus", package = "brassicaData")
    raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:10)])
    raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:10)])
    dat <- intens_theta(raw_napus)
    dat <- remove_suffix(dat, "_Grn")
    dat <- geno_baf_rratio(dat, delthresh = 11)
    expect_equal(class(dat), "norm_data")
    expect_true(is.list(dat))
    expect_true(is.matrix(dat$baf))
    expect_true(is.matrix(dat$rratio))
    expect_true(is.matrix(dat$geno))
    expect_equal(dim(dat$baf), dim(dat$rratio))
    expect_equal(dim(dat$baf), dim(dat$geno))
    expect_equal(dim(dat$baf), c(length(dat$snps), length(dat$samples)))
  }
})