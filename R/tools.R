#' Rename sample names
#'
#' @param dat List object, containing at least two matrices "intensity" and "theta". Or matrix with raw data.
#' @param ... Placeholder for generic function.
#' @return List with two matrices "intensity" (signal intensities) and "theta" (genotype value).
#' @examples
#' if(require(brassicaData)){
#' data(raw_napus, package = "brassicaData")
#' samples <- read_sample_sheets(files = list.files(system.file("extdata",
#' package = "brassicaData"), full.names = TRUE, pattern = "csv"))
#' raw_napus <- rename_samples(raw_napus, samples = samples[,2:1], suffix = c("_Grn", "_Red"))
#' }
#' @export
rename_samples <- function(dat, ...) {
  UseMethod("rename_samples")
}

#' @param samples Matrix or data.frame with two columns: 1. (prefix of) colnames of dat and 2. meaningful sample name.
#' @param suffix Vector of two characters (e.g. c("Grn.idat", "Red.idat")) describing which suffix should be appended to the sample name.
#' @rdname rename_samples
#' @export
rename_samples.raw_data <-
  function(dat, samples, suffix = NULL, ...) {
    samples1 <- samples[, 1]
    samples2 <- samples[, 2]
    cnames <- dat$samples
    if (length(samples[, 1]) == ncol(dat$raw)) {
      #if one name per column
      cnames <- samples2[match(samples1, cnames)]
    }else{
      # if only one name for each pair of columns
      for (i in 1:length(samples1)) {
        coln <- grep(samples1[i], cnames)
        if (!is.null(coln)) {
          if (is.null(suffix)) {
            suffix <- unlist(strsplit(cnames[coln], split = samples1[i]))
            suffix <- suffix[c(2, 4)]
          }
          newName <- paste(samples2[i], suffix , sep = "")
          cnames[coln] <- newName
        }
      }
    }
    dat$samples <- make.names(cnames, unique = TRUE)
    dat
  }

#' @rdname rename_samples
#' @export
rename_samples.norm_data <- function(dat, samples, ...) {
  samples1 <- samples[, 1]
  samples2 <- samples[, 2]
  dat$samples <-
    make.names(samples2[match(samples1, dat$samples)], unique = TRUE)
  dat
}

#' Remove suffix from sample names
#'
#' @param dat List object, containing at least two matrices "intensity" and "theta".
#' @param suffix Vector of  character (e.g. "_Grn" or "Red").
#'
#' @return List with two matrices "intensity" (signal intensities) and "theta" (genotype value).
#' @examples
#' if(require(brassicaData)){
#' data(raw_napus, package = "brassicaData")
#' \dontshow{
#' raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:10)])
#' raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:10)])
#' }
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' }
#' @export
remove_suffix <- function(dat, suffix) {
  dat$samples <- unlist(strsplit(dat$samples,suffix))
  dat
}



#' Find peaks
#'
#' @param dat Vector, containing theta values for one sample and chromosome.
#' @param npeaks Integer, Number of peaks to be detected.
#' @param breaks Integer, Number of breaks for the histogram.
#' @param check Logical, if TRUE it is checked if the central peak is approximately in the middle of the two other peaks.
#' If the data has not heterozygotes, the third peak might be very close to one of the homozygous peaks.
#' In that case the peak is set to the middle of the other two peaks.
#' @param method Character, "mixture", "density" or "histogram".
#' Defines which method is used to transform the data into a distribution.
#' @keywords internal
find_peak <-
  function(dat, npeaks = NULL, breaks = round(length(dat) * 0.1),
           check = FALSE, method = "density") {
    method <-
      match.arg(arg = method, choices = c("mixture", "density", "histogram"))
    if (method == "mixture") {
      require_package("mixtools")
      if (.Platform$OS.type == "unix") {
        sink('/dev/null')
      } else {
        sink("NUL")
      }
      nm2 <-
        suppressMessages(mixtools::normalmixEM(dat, mu = c(min(dat), max(dat))))
      nm3 <-
        suppressMessages(mixtools::normalmixEM(dat, mu = c(
          min(dat), stats::median(dat), max(dat)
        )))
      sink()
      if (nm2$loglik > nm3$loglik) {
        return(nm2$mu)
      }else{
        return(nm3$mu)
      }
    }else if (method == "histogram") {
      histogram <- graphics::hist(dat, breaks = breaks, plot = FALSE)
      y <- histogram$density
      x <- histogram$counts
      z <- histogram$mids
    }else{
      dens <- stats::density(dat)
      y <- dens$y
      x <- dens$y
      z <- dens$x
    }
    peaks <- which(diff(sign(diff(y))) == -2) + 1
    if (is.null(npeaks) | length(peaks) <= npeaks) {
      index <- 1:length(peaks)
    }else if (length(peaks) > npeaks) {
      index <- length(peaks) - npeaks
      index <- sort(x[peaks])[index]
      index <- which(x[peaks] > index)
    }
    peaks <- sort(z[peaks][index])
    if (!is.null(npeaks) & check & length(peaks) == 3) {
      if (npeaks == 3) {
        tmp <- peaks
        tmp <- tmp - tmp[1]
        tmp <- tmp / tmp[3]
        if (tmp[2] < 0.4 | tmp[2] > 0.6) {
          peaks[2] <- (peaks[1] + peaks[3]) / 2
        }
      }
    }
    peaks
  }
#' Require package wrapper
#'
#' Wrapper to require package and return an error, if it is missing.
#' 
#' @importFrom openxlsx read.xlsx
#' @references http://r-pkgs.had.co.nz/description.html#dependencies
#' @keywords internal
require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Please install package '", pkg, "'", call. = FALSE)
  }
}

#' Interpolate cluster means
#'
#' @param theta Vector of theta values.
#' @param hcenters Vector of length three with the horizontal cluster centers.
#' @param vcenters Vector of length three with the vertical cluster centers.
#'
#' @return Numeric of interpolated intensity value.
#'
#' @keywords internal
interpol <- function(theta, hcenters, vcenters) {
  d1 <- abs(theta - hcenters[1])
  d2 <- abs(theta - hcenters[2])
  d3 <- abs(theta - hcenters[3])
  ratio1 <-
    d2 / (d1 + d2) * vcenters[1] + d1 / (d1 + d2) * vcenters[2]
  ratio2 <-
    d2 / (d2 + d3) * vcenters[3] + d3 / (d2 + d3) * vcenters[2]
  ratio <- theta
  ratio[theta < hcenters[2]] <- ratio1[theta < hcenters[2]]
  ratio[theta >= hcenters[2]] <- ratio2[theta >= hcenters[2]]
  ratio[is.null(ratio)] <- mean(vcenters, na.rm = TRUE)
  ratio
}

