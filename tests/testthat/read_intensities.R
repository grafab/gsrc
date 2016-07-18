context("Read raw intensity values")

test_that("read_intensities returns correct output", {
  if(require(brassicaData)){
    files <- list.files(system.file("extdata", package = "brassicaData"),
    full.names = TRUE, pattern = "idat")
    samples <- read_sample_sheets(files = list.files(system.file("extdata",
    package = "brassicaData"), full.names = TRUE, pattern = "csv"),cols = c("Sample_ID", "SentrixBarcode_A", "SentrixPosition_A"))
    column_names <- sapply(strsplit(files,split="/"), FUN=function(x) x[length(x)])
    data("dictionary", package = "brassicaData")
    data("chrPos", package = "brassicaData")
    dat <- read_intensities(files = files, dict = dictionary,
    cnames = column_names, pos = chrPos)
    expect_equal(class(dat), "raw_data")
    expect_true(is.list(dat))
    expect_true(is.vector(dat$chr))
    expect_true(is.vector(dat$pos))
    expect_true(is.vector(dat$snps))
    expect_true(is.vector(dat$samples))
    expect_true(is.matrix(dat$raw))
  }
})