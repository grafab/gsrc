context("Segmentation")

test_that("segm returns correct output", {
  if(require(brassicaData) && require(DNAcopy)){
    data("raw_napus", package = "brassicaData")
    raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:200)])
    raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:20)])
    dat <- intens_theta(raw_napus)
    dat <- remove_suffix(dat, "_Grn")
    dat <- geno_baf_rratio(dat, delthresh = 11)
    dat <- filt_snps(dat, dat$snps[is.na(rowMeans(dat$baf, na.rm = TRUE))])
    set.seed(31415)
    dat <- segm(dat)
    expect_equal(class(dat), "norm_data")
    expect_true(is.list(dat))
    expect_true(is.data.frame(dat$cna))
    expect_equal(dim(dat$cna), c(104,6))
  }
})