## ---- results="hide"-----------------------------------------------------
#library(gsrc)
#require(devtools)
#devtools::install_github("grafab/brassicaData")
#devtools::install_github("grafab/gsrc")
devtools::load_all()

## ---- eval = FALSE-------------------------------------------------------
#  files <- list.files("/YOUR/DATA/REPOSITORY/",
#                      pattern = "idat",full.names = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  files <- list.files(system.file("extdata",
#                                  package = "brassicaData"),
#                      full.names = TRUE,
#                      pattern = "idat")

## ---- eval = FALSE-------------------------------------------------------
#  samples <- read_sample_sheets(files =
#                                  list.files(system.file("extdata",package = "brassicaData"),
#                                             full.names = TRUE,
#                                             pattern = "csv"))

## ---- eval = FALSE-------------------------------------------------------
#  controls <- grep("H2O", samples$Names)
#  if(length(controls) > 0) samples <- samples[-controls, ]
#  files <- grep(paste(samples$ID, collapse = "|"), files, value = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  column_names <- sapply(strsplit(files, split = "/"), FUN=function(x) x[length(x)])

## ------------------------------------------------------------------------
data(dictionary, package = "brassicaData", envir = environment())
head(dictionary)
data(chrPos, package = "brassicaData", envir = environment())
head(chrPos)

## ---- eval = FALSE-------------------------------------------------------
#  raw_data <- read_intensities(files = files,
#                               dict = dictionary,
#                               cnames = column_names,
#                               pos = chrPos)

## ---- eval = FALSE-------------------------------------------------------
#  str(raw_data)

## ---- eval = FALSE-------------------------------------------------------
#  raw_data <- rename_samples(raw_data,
#                             samples = samples[,2:1],
#                             suffix = c("_Grn", "_Red"))

## ------------------------------------------------------------------------
data(raw_napus, package = "brassicaData", envir = environment())

## ---- fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap = "Raw Data Histogram"----
check_raw(raw_napus, thresh = 28000, breaks = 20)

## ---- eval = TRUE--------------------------------------------------------
length(raw_napus$samples)
raw_napus <- filt_samp(raw_napus, check_raw(raw = raw_napus, plot = FALSE, thresh = 28000))
length(raw_napus$samples)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap = "Boxplot comparing green and red signal distibutions"----
boxplot(as.vector(raw_napus$raw[, seq(1, length(raw_napus$samples), 2)]),
        as.vector(raw_napus$raw[, seq(2, length(raw_napus$samples), 2)]),
        names = c("Green", "Red"))

